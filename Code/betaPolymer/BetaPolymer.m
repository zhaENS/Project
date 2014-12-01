classdef BetaPolymer<handle
    % A class for the beta model (beta~=2)
    properties
        params
        betaMatrix
        position
        noise
        beadDist
        handles
        step = 1;        
        savedPosition
        simulationTime 
        encounterHist
        msd
    end
    
    properties (Access=private)
        noiseTerm % saved noise vector, to speed up the system, obj.noise draws values from it 
    end
    
    methods
        
        function obj = BetaPolymer(params)
%             obj.SetDefaultParams;
            if exist('params','var')
                f = fieldnames(params);
                for fIdx = 1:numel(f)
                    if strcmpi(f{fIdx},'beta')
                        if numel(params.beta)==1
                            % set all beta values to the inserted value 
                            params.beta = params.beta*ones(1,params.numBeads);
                        elseif numel(params.beta)~=params.numBeads
                            error('number of beta values should match the number of beads')
                        end
                    end
                    obj.params.(f{fIdx}) = params.(f{fIdx});
                end
            else 
                obj.SetDefaultParams;
            end
            % sort parameters according to the new inserted values
%             obj.params.noiseSTD        = sqrt(2*obj.params.diffusionConst*obj.params.dt);
            
            if (~isempty(obj.params.affineBeadsNum) && strcmpi(obj.params.saveBeadDist,'last'))
                warning('params.saveBeadDist is set to current')
                obj.params.saveBeadDist = 'current';
            end
        end
        
        function SetDefaultParams(obj)
            obj.params.numBeads        = 10;
            obj.params.beta            = 2*ones(1,obj.params.numBeads); % beta values [beta=2=Rouse]
            obj.params.b               = 0.1;
            obj.params.diffusionConst  = 1;
            obj.params.dt              = 1e-5;
            obj.params.dimension       = 3;
            obj.params.noiseSTD        = sqrt(2*obj.params.diffusionConst*obj.params.dt);
            obj.params.numSteps        = 1000;
            obj.params.affineBeadsNum  = []; % two column vector of beads with affinity to one another
            obj.params.kOff            = 0.05; % [some units]
            obj.params.plot            = false;
            obj.params.encounterDist   = obj.params.b/2;
            obj.params.stiffConnectors = [];
            obj.params.connectedBeads  = []; %the beads connected other than the trivial connections. given by bead pairs
            obj.params.recordPath      = false;% save bead positions throughout the simulation (used for later display)
            obj.params.noiseCycle      = 10000; % noise terms will be generated every 1000 steps
            obj.params.saveBeadDist    = 'last';
            obj.params.springConst     = -obj.params.dimension*obj.params.diffusionConst*obj.params.dt/obj.params.b^2;
        end
        
        function Initialize(obj)
            obj.step = 1;
            obj.CreateBetaMatrix
            obj.GetInitialChainPosition;
%             obj.CreateNoise;
        end
        
        function CreateBetaMatrix(obj)
            warning off
            % construct the beta matrix 
            B = zeros(obj.params.numBeads);
            p = 1:obj.params.numBeads-1;
            for lIdx = 1:obj.params.numBeads
                for mIdx = 1:obj.params.numBeads
                    B(mIdx,lIdx) = (4*obj.params.springConst*2/obj.params.numBeads)*...
                        sum((sin(p*pi./(2*obj.params.numBeads)).^obj.params.beta(lIdx)).*cos((lIdx-0.5)*p*pi./...
                        obj.params.numBeads).*cos((mIdx-0.5)*p*pi./obj.params.numBeads));
                end
            end
            
            % Construct the connectivity graph for the polymer 
            gc = diag(ones(1,obj.params.numBeads-1),1)+diag(ones(1,obj.params.numBeads-1),-1);
             for bIdx= 1:size(obj.params.connectedBeads)
                 gc(obj.params.connectedBeads(bIdx,1), obj.params.connectedBeads(bIdx,2))=1;
                 gc(obj.params.connectedBeads(bIdx,2), obj.params.connectedBeads(bIdx,1))=1;
             end
            graphConnectivity = sparse(gc);
            
%             [a,b]             = find(gc);
%             graphConnectivity = sparse(a,b,ones(size(a),1),size(gc,1),size(gc,2));    
            
            % Find all shortest paths in the connectivity graph 
            [sPaths] = graphallshortestpaths(graphConnectivity);            
            sPaths   = sPaths+1; % increas by one to indicate indices
            bTemp    = B;
              % Assign weights according to distance 
            for rIdx = 1:obj.params.numBeads
                if (obj.params.numBeads-rIdx)>= rIdx
                    vals = bTemp(rIdx,rIdx:obj.params.numBeads);
                else
                    vals = bTemp(rIdx,rIdx:-1:1);
                end
                B(rIdx,:) = vals(sPaths(rIdx,:));               
            end
            
            % Change the diagonal elements such that the sum of each row is
            % zero                              
            B = B-diag(diag(B));
            B = B-diag(sum(B,2));
         
            clear bTemp
            obj.betaMatrix = B;
        end
        
        function GetInitialChainPosition(obj)
            obj.position.prev = zeros(obj.params.numBeads,obj.params.dimension);
            r = randn(obj.params.numBeads,3);
            for bIdx = 2:obj.params.numBeads
                obj.position.prev(bIdx,:) = obj.position.prev(bIdx-1,:)+ r(bIdx,:)*0.7*obj.params.b./norm(r(bIdx,:));
            end
%             obj.CreateLoopsByBrownianBridge;
            
            obj.savedPosition(:,:,obj.step) = obj.position.prev;
            obj.GetBeadsDist;
        end
                            
        function Run(obj)
                     
            while obj.step<=obj.params.numSteps
                % advance one step 
              obj.Step
            end            
            obj.CalculateEncounters
        end
        
        function Step(obj)
                % advance one step 
                tic
%                 obj.GetNoise;
                obj.GetBeadsDist;
%                 obj.LinkCloseBeads
                obj.position.cur = obj.betaMatrix*obj.position.prev+...
                                     obj.position.prev+...                                     
                                     randn(obj.params.numBeads,obj.params.dimension)*obj.params.noiseSTD;
%                                  obj.noise;
                ms = [obj.position.cur(:,1)-obj.savedPosition(:,1,1),...
                      obj.position.cur(:,2)-obj.savedPosition(:,2,1),...
                      obj.position.cur(:,3)-obj.savedPosition(:,3,1)];
                obj.msd = sum((ms).^2,2);% msd for present location 
                
                obj.position.prev = obj.position.cur;
                obj.step = obj.step+1;
                if  obj.params.recordPath 
                   obj.savedPosition(:,:,obj.step) = obj.position.prev;
                end
                obj.simulationTime = toc;
        end
        
        function GetBeadsDist(obj)
            % Calculate the pairwise distance between beads
            % This function requires the mex file pdist2mex in the working path
            if strcmpi(obj.params.saveBeadDist,'current')
                obj.beadDist(:,:,1) = pdist2mex([obj.position.prev]',...
                [obj.position.prev]','euc',[],[],[]);
            elseif strcmpi(obj.params.saveBeadDist,'last')
                % keep the first if we want to also perform linking
                if obj.step == obj.params.numSteps || obj.step == 1
                    obj.beadDist(:,:,1) = pdist2mex([obj.position.prev]',...
                    [obj.position.prev]','euc',[],[],[]);
                end
            elseif strcmpi(obj.params.saveBeadDist,'all')                
               obj.beadDist(:,:,obj.step) = pdist2mex([obj.position.prev]',...
                [obj.position.prev]','euc',[],[],[]);
            end
        end
        
        function LinkCloseBeads(obj)% should be revised according to the connectivity 
            % when affine beads are defined, they are connected if in some
            % simulation step the affine beads come closer than the
            % encounter distance.
            if ~isempty(obj.params.affineBeadsNum)
            r = rand(size(obj.params.affineBeadsNum,1),1);
            
            if strcmpi(obj.params.saveBeadDist,'current')
                ind = 1;
            elseif strcmpi(obj.params.saveBeadDist,'all')
                ind = obj.step;
            elseif strcmpi(obj.params.saveBeadDist,'last')
                ind = 1;
            end
            
            for lIdx = 1:size(obj.params.affineBeadsNum,1)
                if obj.beadDist(obj.params.affineBeadsNum(lIdx,1),obj.params.affineBeadsNum(lIdx,2),ind)<obj.params.encounterDist
                    if r(lIdx)>obj.params.kOff*obj.params.dt
                        % change the Rouse matrix to keep beads close
                        obj.betaMatrix(obj.params.affineBeadsNum(lIdx,1),obj.params.affineBeadsNum(lIdx,2))= -1;
                        obj.betaMatrix(obj.params.affineBeadsNum(lIdx,2),obj.params.affineBeadsNum(lIdx,1))= -1;
                        % adjust diagonal entries
                        obj.betaMatrix(1:obj.params.numBeads+1:obj.params.numBeads^2) = sum(obj.betaMatrix==-1);
                    else
                        % diconnect the beads
                        obj.betaMatrix(obj.params.affineBeadsNum(lIdx,1),obj.params.affineBeadsNum(lIdx,2))= 0;
                        obj.betaMatrix(obj.params.affineBeadsNum(lIdx,2),obj.params.affineBeadsNum(lIdx,1))= 0;
                        % adjust diagonal entries
                        obj.betaMatrix(1:obj.params.numBeads+1:obj.params.numBeads^2) = sum(obj.betaMatrix==-1);
                    end
                end
            end
            end
        end
        
        function Reset(obj)
            obj.step = 1;
        end
        
        function CreateNoise(obj)
            % create noise terms at the begining of simulation to save
            % running time 
            m = min([obj.params.noiseCycle, obj.params.numSteps]);
            
            obj.noiseTerm = randn(obj.params.numBeads,obj.params.dimension,m)*obj.params.noiseSTD;
        end
        
        function GetNoise(obj)
            noiseInd = mod(obj.step,obj.params.noiseCycle);
            if noiseInd==0
                obj.CreateNoise
                noiseInd = obj.params.noiseCycle;
            end
            obj.noise = obj.noiseTerm(:,:,noiseInd);
            obj.noise(obj.params.stiffConnectors,:) = 0; 
        end
        
        function CalculateEncounters(obj)
            obj.encounterHist = sum(obj.beadDist<obj.params.encounterDist,3);
            obj.encounterHist = obj.encounterHist-diag(ones(1,obj.params.numBeads)).*obj.encounterHist;
%             for bIdx = 1:numel(obj.params.stiffConnectors)
%             obj.encounterHist
%             end
        end
    end
end