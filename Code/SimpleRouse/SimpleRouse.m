classdef SimpleRouse<handle
    
    properties
        params
        rouseMatrix
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
        
        function obj = SimpleRouse(params)
%             obj.SetDefaultParams;
            if exist('params','var')
                f = fieldnames(params);
                for fIdx = 1:numel(f)
                    obj.params.(f{fIdx}) = params.(f{fIdx});
                end
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
            obj.params.b               = 0.1;
            obj.params.diffusionConst  = 1;
            obj.params.dt              = 1e-5;
            obj.params.dimension       = 3;
            obj.params.noiseSTD        = sqrt(2*obj.params.diffusionConst*obj.params.dt);
            obj.params.numSteps        = 50000;
            obj.params.affineBeadsNum  = []; % two column vector of beads with affinity to one another
            obj.params.kOff            = 0.05; % [some units]
            obj.params.plot            = false;
            obj.params.encounterDist   = obj.params.b/2;
            obj.params.stiffConnectors = [];
            obj.params.connectedBeads  = []; %the beads connected other than the trivial connections. given by bead pairs
            obj.params.recordPath      = false;% save bead positions throughout the simulation (used for later display)
            obj.params.noiseCycle      = 10000; % noise terms will be generated every 1000 steps
            obj.params.saveBeadDist    = 'last';
        end
        
        function Initialize(obj)
            obj.step = 1;
            obj.CreateRouseMatrix
            obj.GetInitialChainPosition;
%             obj.CreateControls
%             obj.CreateNoise;
        end
        
        function CreateRouseMatrix(obj)
            R = zeros(obj.params.numBeads);
            R = R+ diag(-1*ones(1,obj.params.numBeads-1),-1) +...
                diag(-1*ones(1,obj.params.numBeads-1),1);
            for bIdx = 1:size(obj.params.connectedBeads,1)
                R(obj.params.connectedBeads(bIdx,1),obj.params.connectedBeads(bIdx,2)) = -1;
                R(obj.params.connectedBeads(bIdx,2),obj.params.connectedBeads(bIdx,1)) = -1;
            end
            R(obj.params.stiffConnectors,:) = 0;
            R = R+diag(sum(R==-1,2));
            
% %             if half of the matrix are zeroes, switch to sparse
% % %             representation
%             if obj.params.numBeads>20;% numel(find(R==0))/obj.params.numBeads^2>0.5
%               obj.rouseMatrix = sparse(obj.params.springConst*R);
%             else
              obj.rouseMatrix = obj.params.springConst*(R);
%             end
        end
        
        function GetInitialChainPosition(obj)
            obj.position.prev = zeros(obj.params.numBeads,obj.params.dimension);
            r = randn(obj.params.numBeads,3);
            for bIdx = 2:obj.params.numBeads
                obj.position.prev(bIdx,:) = obj.position.prev(bIdx-1,:)+ r(bIdx,:)*0.7*obj.params.b./norm(r(bIdx,:));
            end
            obj.msd = zeros(obj.params.numBeads,1);
%             obj.CreateLoopsByBrownianBridge;
            
            obj.savedPosition(:,:,obj.step) = obj.position.prev;
            obj.GetBeadsDist;
        end
        
        function CreateLoopsByBrownianBridge(obj)%[unfinished]
            % if connected beads are define, change bead locations by
            % bridging the connected beads by a Brownian beads. This should
            % save some running time 
            if ~isempty(obj.params.connectedBeads)
                % sort the loops according to size
                % start with the biggest loop and sequentially bridge all
                % connectors 
                loopSize = obj.params.connectedBeads(:,2)-obj.params.connectedBeads(:,1);
                [~, loopInds] = sort(loopSize,'Descend');
                for lIdx = 1:numel(loopsInds)
                    posStart = obj.position.prev(obj.params.connectedBeads(loopInds(lIdx),1),:);
                    posEnd   = obj.position.prev(obj.params.connectedBeads(loopInds(lIdx),2),:);
                    % replace all the path between the connected beads 
 
                end
                
            end
        end
            
        function CreateControls(obj)% should be removed (obsolete)
            if obj.params.plot
                obj.handles.graphical.mainFigure = figure('Name','RouseSimulation',...
                    'Units','norm');
                xMin = min(obj.position.prev(:,1));
                xMax = max(obj.position.prev(:,1));
                yMin = min(obj.position.prev(:,2));
                yMax = max(obj.position.prev(:,2));
                zMin = min(obj.position.prev(:,3));
                zMax = max(obj.position.prev(:,3));
                k = 10*obj.params.b;
                obj.handles.graphical.mainAxes = axes('Parent',obj.handles.graphical.mainFigure,...
                    'Units','norm',...
                    'XLim',[xMin-k xMax+k],...
                    'YLim',[yMin-k yMax+k],...
                    'ZLim',[zMin-k, zMax+k]);
                % plot initial chain positino
                obj.handles.graphical.chain = line('XData',obj.position.prev(:,1),...
                    'YData',obj.position.prev(:,2),...
                    'ZData',obj.position.prev(:,3),...
                    'Parent',obj.handles.graphical.mainAxes);
                cameratoolbar(obj.handles.graphical.mainFigure);
            end
        end
        
        function Run(obj)
                        
            while obj.step<=obj.params.numSteps
              obj.Step
            end
            
            obj.CalculateEncounters
        end
        
        function Step(obj)
                % advance one step 
                 tic;
%                 springConst = -(obj.params.dimension*obj.params.diffusionConst*obj.params.dt/(obj.params.b^2));
%                 obj.GetNoise;
                obj.GetBeadsDist;
%                 obj.LinkCloseBeads
                obj.position.cur = obj.rouseMatrix*obj.position.prev+...
                                     obj.position.prev+...                                     
                                     randn(obj.params.numBeads,obj.params.dimension)*obj.params.noiseSTD;
                ms = [obj.position.cur(:,1)-obj.savedPosition(:,1,1),...
                      obj.position.cur(:,2)-obj.savedPosition(:,2,1),...
                      obj.position.cur(:,3)-obj.savedPosition(:,3,1)];
                obj.msd = sum((ms).^2,2);% msd for present location 
                obj.position.prev = obj.position.cur;
                
                obj.step = obj.step+1;
                if  obj.params.recordPath 
                   obj.savedPosition(:,:,obj.step) = obj.position.prev;
                end
%                 t2  = clock; 
                obj.simulationTime = toc;% etime(t2,t1);
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
        
        function LinkCloseBeads(obj)
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
                        obj.rouseMatrix(obj.params.affineBeadsNum(lIdx,1),obj.params.affineBeadsNum(lIdx,2))= -1;
                        obj.rouseMatrix(obj.params.affineBeadsNum(lIdx,2),obj.params.affineBeadsNum(lIdx,1))= -1;
                        % adjust diagonal entries
                        obj.rouseMatrix(1:obj.params.numBeads+1:obj.params.numBeads^2) = sum(obj.rouseMatrix==-1);
                    else
                        % diconnect the beads
                        obj.rouseMatrix(obj.params.affineBeadsNum(lIdx,1),obj.params.affineBeadsNum(lIdx,2))= 0;
                        obj.rouseMatrix(obj.params.affineBeadsNum(lIdx,2),obj.params.affineBeadsNum(lIdx,1))= 0;
                        % adjust diagonal entries
                        obj.rouseMatrix(1:obj.params.numBeads+1:obj.params.numBeads^2) = sum(obj.rouseMatrix==-1);
                    end
                end
            end
            end
        end
        
        function Plot(obj)
            if obj.params.plot
                set(obj.handles.graphical.chain,'XData',obj.position.prev(:,1),...
                    'YData',obj.position.prev(:,2),...
                    'ZData',obj.position.prev(:,3));
                drawnow
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