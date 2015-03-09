classdef CalculateBeadDistancesByRouseModel<handle
    
    properties (Access=public)
        encounterMat
        connectivityMat
        params           % class for parameters
        graph
        chain           % the Rouse chain class
    end
    
    properties (Access=private)
        smoother =Smoother; % signal smoother class
    end
    
    methods
        
        function obj = CalculateBeadDistancesByRouseModel()
            % class constructor 
         obj.params = ReconstructionParams;
%             obj.SetParams(params);
        end
                
        function Initialize(obj,encounterMat)
            % The input should be two sided encounter ~histogram~ matrix from the HiC
            % experiment 
            % two sided analysis should then follow
            % the connectivity graph should  be created based on two sided
            % analysis 
             if size(encounterMat,1)~=((size(encounterMat,2)-1)/2+1)
                 error('encounter matrix must reflect two sides of a square matrix and must be of size Nx(2*N+1)')
             end
             
             obj.ProcessEncounterMatrix(encounterMat);
            
             obj.AnalyzeConnectivity

            
            %  create a chain from the graph
            obj.CreateChainFromConnectivityGraph;
            obj.DisplayChain;
        end       
        
        function AnalyzeConnectivity(obj)
             % Perform analysis for each side        
             eMatLeft  = obj.TransformEncountersToConnectivityMap('left');
             eMatRight = obj.TransformEncountersToConnectivityMap('right');                                      
             eMat      = eMatLeft |eMatRight;
             obj.connectivityMat = eMat;
%             obj.DisplayConnectivityGraph(eMat,above,obj.distToAnalyze,obj.beadsToAnalyze);
        end
        
        function eMat= TransformEncountersToConnectivityMap(obj,side)
            % Analyze the encounters and determine bead neighbors in the
            % specified distance
            if strcmpi(side,'left')% left encounters
                % take the left side of the encounter matrix and rotate it 
%                 encounters = obj.encounterMat(obj.params.reconstruction.beadRange.bead1,...
%                                               1:(obj.params.reconstruction.beadRange.bead2(end)-1)/2);
%                 encounters = fliplr(encounters);
                 encounters = obj.encounterMat(:,1:(size(obj.encounterMat,2)-1)/2);% tril(obj.encounterMat)';
                 encounters = fliplr(encounters);
            elseif strcmpi(side,'right')% right encounters 
%                 encounters = obj.encounterMat(obj.params.reconstruction.beadRange.bead1,(obj.params.reconstruction.beadRange.bead2(end)+1)/2 +1 :end);                
                  encounters = obj.encounterMat(:,(size(obj.encounterMat,2)+1)/2:end);
            end
            
            % Preallocations
            above = cell(1,size(encounters,1));% save indices of distances falling above the nearest neighor encounter probability
            dists = cell(size(encounters,1),size(encounters,2));
            histK = cell(size(encounters,1),size(encounters,2));
            chainRingStruct = cell(size(encounters,1),1);
            % Construct a binary connection matrix for a specific distance
            eMat = false(size(encounters,1),...
                         size(encounters,1),...
                         obj.params.reconstruction.numDistances);
            di          = diag(ones(1,size(eMat,1)-1),1)+diag(ones(1,size(eMat,1)-1),-1);% include nearest neighbors by default
            chainRingImage = zeros(size(encounters,1), size(encounters,2)*2 +1,3);  % an image representing positions of chains and rings                        
            for bIdx = 1:size(encounters,1);% for each bead
                
                  if strcmpi(side,'left')              
                  observedProb = encounters(bIdx,1:(bIdx-1));
                  else
                      observedProb = encounters(bIdx,1:(size(encounters,2)-bIdx+1));
                  end
                  
                if ~all(isnan(observedProb))                                        
                    % Divide the probabilites into distances according to a division given by the
                    % expected model

                    [~,modelValues]         = obj.TransformProbToDist(observedProb);
%                     [chainRingStruct{bIdx}] = obj.AnalyzeEncounterAsRingsAndChains(observedProb);
                    % fill in the ring chain image 
                    for cIdx = 1:numel(chainRingStruct{bIdx});
                        if strcmpi(chainRingStruct{bIdx}(cIdx).type,'chain')
                            if strcmpi(side,'left')
                              chainRingImage(bIdx,(chainRingStruct{bIdx}(cIdx).startInd:chainRingStruct{bIdx}(cIdx).endInd),1)=255;
                            elseif strcmpi(side,'right')
                              chainRingImage(bIdx,chainRingStruct{bIdx}(cIdx).startInd:chainRingStruct{bIdx}(cIdx).endInd,1)=255;  
                            end
                        elseif strcmpi(chainRingStruct{bIdx}(cIdx).type,'ring')
                            if strcmpi(side,'left')
                              chainRingImage(bIdx,(chainRingStruct{bIdx}(cIdx).startInd:chainRingStruct{bIdx}(cIdx).endInd),3)=255;
                            elseif strcmpi(side,'right')
                              chainRingImage(bIdx,chainRingStruct{bIdx}(cIdx).startInd:chainRingStruct{bIdx}(cIdx).endInd,3)=255;  
                            end
                        end
                    end
                    % obtain the values of the composite structure
%                     vals = obj.GetCompositeFunctionVals(chainRingStruct);
                    inds = 1:numel(observedProb);
                    % display
                    if mod(bIdx,400)==0
                        obj.PlotBeadClusteringByDistance(observedProb,inds,modelValues);
                        title(num2str(bIdx))
                    end
                                                          
                    % Calculate the histogram
                    above{bIdx} = find(observedProb>modelValues(1));
%                     below{bIdx} = [];
                    for kIdx = 1:numel(modelValues)-1
                        dists{bIdx,kIdx} = find(observedProb>modelValues(kIdx+1) & observedProb<=modelValues(kIdx));
%                         below{bIdx,kIdx} = dists{bIdx,kIdx}((dists{bIdx,kIdx}<kIdx));                                                
                        % add the terms "above" to the dist 1 neighbors
                        histK{bIdx,kIdx} = numel(dists{bIdx,kIdx});
                        if kIdx ==1
                            dists{bIdx,kIdx} = [dists{bIdx,kIdx}, above{bIdx}];
                        end
                    end                    
                end
            end
            
            for dIdx = obj.params.reconstruction.distToAnalyze
                eMat(:,:,dIdx) = eMat(:,:,dIdx)|di;
                for b1Idx = 1:size(encounters,1)
                    % Collect all beads at distance 1
                    inds1 = b1Idx+dists{b1Idx,dIdx};
                    inds1 = inds1(inds1<size(encounters,2));
                    inds2 = b1Idx-dists{b1Idx,dIdx};
                    inds2 = inds2(inds2>=1);
                    if strcmpi(side,'left')
                    eMat(b1Idx,inds2,dIdx)= true;
                    elseif strcmpi(side,'right')
                        eMat(b1Idx,inds1,dIdx)= true;
                    end
                end
            end
            
        end
                
        function [expectedSignal,fitStructModel] = GetExpectedEncounterSignal(obj,encounterMat)
            % Get the expected encounter signal from an encounter matrix 
            expectedSignal    = obj.MeanIgnoreNaN(encounterMat);
            expectedSignal    = obj.smoother.Smooth(expectedSignal',obj.param.smoothing.method,obj.params.smoothing.nHoodRad,obj.params.smoothing.sigma);
            se                = obj.SumIgnoreNaN(expectedSignal(:,1,1));% normalize
            expectedSignal    = expectedSignal(:,1,1)./se;
            [fitStructModel]  = fit((1:numel(expectedSignal))',expectedSignal,obj.params.optimization.model,obj.params.optimization.fitOpt);
        end        
        
        function GetDistanceDistribution(obj,dist)
            % calculate the distribution for a specifind distance
            f = figure;
            a = axes('Parent',f,'FontSize',30);
            for dIdx = 1:numel(dist)
                d     = obj.encounterMat(:,dist(dIdx));
                d     = d(~isnan(d));
                d     = d(d~=0);
                d     = (d-mean(d));
                s     = std(d);
                if s~=0
                    d=d./s;
                    
                    [v,e] = ecdf(d);
                    line('XData',e,...
                        'YData',v,...
                        'DisplayName',num2str(dist(dIdx)),...
                        'Color',rand(1,3),...
                        'LineWidth',4,...
                        'Parent',a);
                end
                
            end
            l=legend(get(a,'Children'));
            set(l,'FontSize',10);
            xlabel(a,'(Distance-\mu)/\sigma');
            ylabel(a,'CDF');
            
        end
        
        function DisplayConnectivityGraph(obj,eMat,above,distToAnalyze,beadToAnalyze)
            % construct a graph
            if ~exist('beadToAnalyze','var')
                beadToAnalyze = 1:size(eMat,1);
            end
            
            obj.connectivityMat         = triu(eMat(:,:,distToAnalyze));
%             inds                        = setdiff(1:size(eMat,1),beadToAnalyze);
%             obj.connectivityMat(inds,:) = false;
            % add nearest neighbor connectivity 
            obj.connectivityMat = obj.connectivityMat | diag(true(1,size(eMat,1)-1),1);
            obj.graph                   = biograph(obj.connectivityMat);
            set(obj.graph,'LayoutType','hierarchical','EdgeType','straight','NodeCallback',@obj.NodeCallback);
            
            % mark edges between nodes that have higher probability than nearest
            % neighbor with red
            for aIdx = 1:size(eMat,1)
                set(obj.graph.Nodes(aIdx),'Label',['Bead ' num2str(aIdx)]);
                for a1Idx = 1:numel(above{aIdx})
                    sourceNode = ['Node ' num2str(aIdx)];
                    if (aIdx +a1Idx)<=numel(obj.beadRange.bead2)
                        sinkNode =  ['Node ' num2str(aIdx+a1Idx)];
                    else
                        sinkNode = ['Node ' num2str(aIdx-a1Idx)];
                    end
                    
                    f = obj.graph.getedgesbynodeid(sourceNode,sinkNode);
                    set(f,'LineColor',[1 0 0]);
                end
            end
            view(obj.graph);
            
        end
        
        function CreateChainFromConnectivityGraph(obj)
            connectedBeads = obj.connectivityMat(:,:,1);% always analyze the first neighbor map 
            % remove the trivial connections on the main diagonal
            connectedBeads = triu(connectedBeads-diag(diag(connectedBeads,1),1)); 
            [r,c]       = find(connectedBeads);    
            obj.params.chain.numBeads =size(obj.encounterMat,1);
            obj.params.chain.connectedBeads = [r,c];                        
            obj.params.chain.springConst     = -(obj.params.chain.dimension* obj.params.chain.diffusionConst*...
                                                 obj.params.chain.dt/obj.params.chain.b^2)*obj.connectivityMat;% ones(obj.params.chain.numBeads); % can be a scalar or a matrix the size of (numBeads) X (numBeads)         
            obj.chain = SimpleRouseFramework;
            obj.chain.Initialize(obj.params.chain);
%             obj.chain.Run;
            
        end
        
        function DisplayChain(obj)
            % display the connected chain 
            ChainDynamicsPlayer(obj.chain);
        end
        
    end
    
    methods (Access=private)
        
        function ProcessEncounterMatrix(obj,encounterMat)
            % Interpolate and normalize the encounter histogram             
            obj.encounterMat                          = encounterMat;
            obj.params.reconstruction.beadRange.bead1 = 1:size(encounterMat,1);
            obj.params.reconstruction.beadRange.bead2 = 1:size(encounterMat,2);
            
            obj.encounterMat(isnan(obj.encounterMat))=0;% zero out nans

            % Smooth left and right parts of the encounter matrix                
            left  = fliplr(obj.encounterMat(:,1:(size(encounterMat,2)-1)/2));
            
            right = (obj.encounterMat(:,(size(encounterMat,2)+1)/2 +1 :end));
            
            regOrder = 1; % regularization order [0,1,2]
            lambda   = 0.01; % regularization constant
            minNumPts = 4; % min number of points to perform analysis (smoothing)
            % calculate the sum for each row to be used for the
            % normalization of both sides

            

            % interpolate zero values in the signal            
            for lIdx = minNumPts:size(obj.encounterMat,1)
                left(lIdx,1:lIdx-1) = obj.InterpolateZeroValuesInSignal(left(lIdx,1:lIdx-1));
                left(lIdx,lIdx:end) = 0; % zero out the ends

            end
%             ml    = obj.MaxIgnoreNaN(left(:));% maxima of the left signal 
            
            right(isnan(right)) = 0;
            for rIdx = 1:size(right,1)-minNumPts
                right(rIdx,1:end-rIdx+1)     = obj.InterpolateZeroValuesInSignal(right(rIdx,1:end-rIdx+1));
                right(rIdx,end-rIdx+2:end)   = 0; 
            end  
            
            s = sum(right+left,2);% normalization constant for each row
            for sIdx = 1:numel(s)
                if s(sIdx)>0
                    left(sIdx,:) = left(sIdx,:)./s(sIdx);
                    right(sIdx,:) = right(sIdx,:)./s(sIdx);
                end            
            end
            % perform smoothing
            for lIdx = minNumPts:size(obj.encounterMat,1)
               if ~all(left(lIdx,:)==0)% perform smoothing for non zero signals
                % smooth using the inverse heat equation
                
                  [~,left(lIdx,1:lIdx-1)] =BoundaryElementHeatEquation('TestBemHeatEq_optimized',left(lIdx,1:lIdx-1),regOrder,lambda,false); 
               end
            end
            
            for rIdx =1:size(right)-minNumPts
               if ~all(right(rIdx,:)==0)
                  [~,right(rIdx,1:end-rIdx+1)] = BoundaryElementHeatEquation('TestBemHeatEq_optimized',right(rIdx,1:end-rIdx+1),regOrder,lambda,false); 
               end
            end
            
%             mr    = obj.MaxIgnoreNaN(right(:));% maxima of the right signal
%             
%             obj.smoother.Smooth(left./ml,obj.params.smoothing.method,obj.params.smoothing.nHoodRad, obj.params.smoothing.sigma,obj.params.smoothing.kernel);
%             left  = obj.smoother.signalOut.*ml;
%             obj.smoother.Smooth(right./mr,obj.params.smoothing.method,obj.params.smoothing.nHoodRad, obj.params.smoothing.sigma,obj.params.smoothing.kernel);
%             right = obj.smoother.signalOut.*mr;
            obj.encounterMat = [fliplr(left), zeros(size(left,1),1), right];
            
            % Normalize the smoothed encounter matrix
            for bIdx=obj.params.reconstruction.beadRange.bead1
%                 obj.encounterMat(bIdx,:) = obj.InterpolateZeroValuesInSignal(obj.encounterMat(bIdx,:)); % interpolate zero values
                obj.encounterMat(bIdx,:)= obj.encounterMat(bIdx,:)./obj.SumIgnoreNaN(obj.encounterMat(bIdx,:)); % normalize to get probabilities 
            end
        end
        
        function [sigOut] = InterpolateZeroValuesInSignal(obj,sigIn)
            % remove nan values by interpolation of the signal
            sigOut = sigIn;
            s      = size(sigIn);
            if any(s(1:2)==1) % for 1D signal
                zeroInds   = find(sigIn==0 |isnan(sigIn));
                noZeroInds = find(~(sigIn==0) & ~isnan(sigIn));
                if ~isempty(zeroInds) && numel(noZeroInds)>2
                    % Interpolate the signal in the nan positions
                    x     = noZeroInds;
                    y     = sigIn(x);
                    sigOut(zeroInds)= interp1(x,y,zeroInds,obj.params.interpolation.zeroInterpolationMethod);% for the boundary values
                    % extrapolate the end values
                    
                    f =find(~isnan(sigOut),1,'first');
                    if f~=1
                        sigOut(1:f) = sigOut(f);
                    end
                    f= find(~isnan(sigOut),1,'last');
                    if f~=numel(sigOut)
                        sigOut(f:end) = sigOut(f);
                    end                                        
                end
            else
                %                 sigOut(isnan(sigIn))=0;
                [zeroInds(:,1), zeroInds(:,2)]     = find(sigOut==0);
                [noZeroInds(:,1), noZeroInds(:,2)] = find(~(sigIn==0) & ~isnan(sigIn));
                sigNoZero = sigIn(:);
                sigNoZero = sigNoZero (sigNoZero ~=0 & ~isnan(sigNoZero));
                if ~isempty(zeroInds)
                    % Interpolate the signal in the nan positions
                    %                     [x,y] = meshgrid(1:size(sigOut,1),1:size(sigOut,2));
                    intPoints = interp2(noZeroInds(:,2), noZeroInds(:,1),sigOut,zeroInds(:,2), zeroInds(:,1),obj.params.interpolation.zeroInterpolationMethod);
                    s         = sub2ind(size(sigOut),zeroInds(:,1), zeroInds(:,2));
                    sigOut(s) = intPoints;
                end
            end
        end
        
        function [chainRingStruct] = AnalyzeEncounterAsRingsAndChains(obj,prob)
            % Seperate the encounter probability signal prob
            % into regions of chains, and rings. The decomposition is
            % defined by the position of the local maximas of the signal
            % prob.
            
            chainStruct  = struct('type',[],...
                                  'containedIn',[],...
                                  'startInd',[],...
                                  'endInd',[],...
                                  'equation',[],...
                                  'length',[],...
                                  'normalizationConst',[],...
                                  'containing',[]);
            ringStruct        = chainStruct;
            minProbVecLength  = 2;
            if numel(prob)>minProbVecLength % analyze only vectors with more than minProbVecLength values
                
            % Find local max in prob signal
            [lMax]       = local_max(prob);
            
            % If index 1 exists as a local_max, remove it
            lMax = lMax(lMax~=1);
            
            % For each max point, find the first point of intersection to
            % its left on the probability signal curve
            % Match position of the NaNs in the signal
                x = 1:numel(prob);
                x(isnan(prob))= nan;
                
            % Start with all rings
            for lmIdx = 1:numel(lMax)
                                
                intersections   = polyxpoly(x,prob,1:numel(prob), prob(lMax(lmIdx)).*ones(1,numel(prob)));
                intersections   = round(intersections); % round to get indices
                % Find the first intersection index to the left of the local max
                d               = find(intersections<lMax(lmIdx),1,'last');
                if isempty(d) && prob(intersections(1))>=max(prob(1:min(numel(prob),10)))                    
                    ringStruct(lmIdx).startInd  = 1;% assume it intersects the first bead
                else
                    ringStruct(lmIdx).startInd  = intersections(d);
                end
                
                ringStruct(lmIdx).endInd    = lMax(lmIdx);% the local max index itself
                ringStruct(lmIdx).length    = lMax(lmIdx)-intersections(d);
                ringStruct(lmIdx).equation  = @(d,N,sig)(sig+(d./N).*(N-d)).^(-1.5); % where sig is the sigma for the containing structure
                ringStruct(lmIdx).type      = 'ring';
%                 ringStruct(lmIdx).normalizationConst = (sum(ringStruct(lmIdx).equation(1:ringStruct(lmIdx).length-1,ringStruct(lmIdx).length,0)));
            end
            
            
            % place chains where there are no rings 
            
            % Sort by ring size
            [loopSize,inds] = sort([ringStruct.length],'descend');
            
            % Rearrange ring structures according to loop size
            ringStruct = ringStruct(inds);
            
            % Construct a matrix with a visual display of the loops
            loopChainMat = zeros(2*numel(ringStruct),numel(prob));
            for lIdx = 1:numel(loopSize);
                loopChainMat(2*lIdx-1,ringStruct(lIdx).startInd:ringStruct(lIdx).endInd)=lIdx;
            end
            
            % Set the chain positions where there are no loops
            l = sum(loopChainMat);
            c = l==0;
            r = regionprops(c,'pixelList');
            
            % Expand the loopChainMat to allow chain insertion (+ 1 line
            % just to make sure there is no overlap)
            loopChainMat = [loopChainMat;2*zeros(numel(r)+1,size(loopChainMat,2))];
            
            for rIdx = 1:numel(r)
                chainStruct(rIdx).startInd = r(rIdx).PixelList(1,1);
                chainStruct(rIdx).endInd   = r(rIdx).PixelList(end,1);
                chainStruct(rIdx).equation = @(d,sig) (sig+d).^(-1.5); %where sig is the sigma of the containing structure
                chainStruct(rIdx).length   = chainStruct(rIdx).endInd-chainStruct(rIdx).startInd;
                chainStruct(rIdx).type     = 'chain';
                chainStruct(rIdx).normalizationConst = sum(chainStruct(rIdx).equation(1:chainStruct(rIdx).length,0));
                loopChainMat(end-2*rIdx+1,chainStruct(rIdx).startInd:chainStruct(rIdx).endInd) = numel(ringStruct)+rIdx;
            end
            
            % Find the correct ordering of loops and chains by labeling
            rp   = regionprops(logical(loopChainMat),'PixelList');
            oMap = zeros(numel(rp),2);
            
            % first column is the order in the composite structure
            % second column is the index in the chainLoopMat
            for rIdx = 1:numel(rp)
                oMap(rIdx,1) = rIdx;
                oMap(rIdx,2) = loopChainMat(rp(rIdx).PixelList(1,2),rp(rIdx).PixelList(1,1));
            end
            % Create the mapping between the ring-chain segments and the
            % order of the composite structures
            chainRingStruct  = [ringStruct,chainStruct];
            chainRingStruct  = chainRingStruct(oMap(:,2));
%             rp = rp(oMap(:,1));
            % Assign the containning order
            lcMat = loopChainMat;
            for oIdx = 1:numel(oMap(:,2))
                lcMat(loopChainMat==oMap(oIdx,2))=oMap(oIdx,1);
            end
            loopChainMat = lcMat;
            for oIdx = 1:numel(oMap(:,2))
                containedIn = find(loopChainMat(1:rp(oIdx).PixelList(1,2)-1,rp(oIdx).PixelList(1,1)),1,'last');
                if ~isempty(containedIn)
                    chainRingStruct(oIdx).containedIn = loopChainMat(containedIn,rp(oIdx).PixelList(1,1));
                end
                
                ind = 1 ;
                for pIdx = 1:numel(rp(oIdx).PixelList(:,1))% go over all pixels of that structure and search downward for the first structure
                 containing = find(loopChainMat(rp(oIdx).PixelList(pIdx,2)+1:end,rp(oIdx).PixelList(pIdx,1)),1,'first');
                 if ~isempty(containing)
                  chainRingStruct(oIdx).containing(ind) = loopChainMat(rp(oIdx).PixelList(pIdx,2)+containing,rp(oIdx).PixelList(pIdx,1));
                  ind = ind+1;
                 end
                end
                 chainRingStruct(oIdx).containing= unique(chainRingStruct(oIdx).containing);
            end
            
            else
                % By default, if the prob vector is shorter than minProbVecLength, set it
                % to be a chain 
                chainRingStruct(1)          = chainStruct;
                chainRingStruct(1).type     = 'chain';
                chainRingStruct(1).startInd = 1;
                chainRingStruct(1).endInd   = numel(prob);
                chainRingStruct(1).length   = numel(prob);
                
            end
        end
        
        function vals = GetCompositeFunctionVals(obj,chainRingStruct)%unfinished
            % Calculate the values of the composite ring-chain structure;
            vals=[];
            for sIdx = 1:numel(chainRingStruct)
                % Evaluate each one of the functions
                
            end
        end
        
        function [dists, modelValues] = TransformProbToDist(obj,prob)
            % Transform the probabilitiy observed into distances
            if strcmpi(obj.params.reconstruction.prob2distMethod,'fitModel')
                % Fit beta to a rouse-like model 
                inds        = find(~isnan(prob));
                [fitStruct] = fit(inds',prob(inds)',obj.params.optimization.model,obj.params.optimization.fitOpt);
                beta        = fitStruct.beta;
                s           = sum(inds.^(-beta));
                modelValues = obj.params.optimization.model(beta,inds);
                modelValues = modelValues./modelValues(1) *max(prob(1:5));
                dists       = (prob *s).^(-1./beta);
            elseif strcmpi(obj.params.reconstruction.prob2distMethod,'rouse')
                % use rouse with beta=1.5
                beta        = 1.5;
                inds        = find(~isnan(prob));
                s           = sum(inds.^(-beta));
                modelValues = obj.params.optimization.model(beta,inds);
                modelValues = modelValues./modelValues(1) *max(prob(1:min(15,size(prob,2))));
                dists       = (prob *s).^(-1/beta);
            elseif strcmpi(obj.params.reconstruction.prob2distMethod,'composite')
%                 do nothing
            end
        end
        
        function NodeCallback(obj,varargin)
            disp('node')
            
            node = varargin{1};
            if strcmpi(node.UserData,'On')
                node.UserData = 'Off';
                nodeColor = [1 1 0.7000];
            else
                node.UserData = 'On';
                nodeColor = [0 1 0];
            end
            a = node.getancestors;
            d = node.getdescendants;
            for aIdx = 1:numel(a);
                a(aIdx).Color = nodeColor;
                a(aIdx).hgUpdate;
            end
            
            for dIdx = 1:numel(d)
                d(dIdx).Color=nodeColor;
                d(dIdx).hgUpdate;
            end
        end
        
        function [observedProb] = ProcessBeadEncounterSignal(obj,encounterSignal)
            % process a single encounter signal
            % fill NaN position with nearest neighbors mean value
            observedProb = obj.InterpolateZeroValuesInSignal(encounterSignal);
            observedProb = obj.smoother.Smooth(observedProb,obj.params.smoothing.method,obj.params.smoothing.nHoodRad,obj.params.smoothing.sigma);            
            sop            = obj.SumIgnoreNaN(observedProb(1,:,end));
            observedProb   = observedProb(1,:,end)./sop; % normalize
        end
    end
    
    methods (Static)
        
        function m = MeanIgnoreNaN(sigIn)
            sigIn(isnan(sigIn)) = 0;
            m = mean(sigIn);
            m = m./sum(m);
        end
        
        function s = SumIgnoreNaN(sigIn)
            % calculate the sum of an input signal without the NaNs
               s  = sum(sigIn(~isnan(sigIn)));
        end
        
        function m = MaxIgnoreNaN(sigIn)
            s = sigIn(:);
            s(isnan(s))= -Inf;
            m = max(s);
        end
        
        function PlotBeadClusteringByDistance(observedProb,inds, k)
            figure,
            plot(inds,observedProb,'bo-',1:numel(k),k,'r','Linewidth',4,'MarkerSize',6),
            set(gca,'FontSize',35,'NextPlot','Add');%,'XScale','log','YScale','log'),
            xlabel('Distance [beads]'),
            ylabel('encoutner Prob.')
            
            % add patches to represent the  distance by probability, given the model
            cMap = rand(numel(inds),3);
            for kIdx = 1:numel(inds)-1
                patch([inds(1) inds(end), inds(end), inds(1)], [k(kIdx) k(kIdx), k(kIdx+1), k(kIdx+1)],...
                    'r','FaceAlpha',0.25,'FaceColor',cMap(kIdx,:));
            end
            set(gca,'XLim',[1 inds(end)])
        end
        
    end
    
end