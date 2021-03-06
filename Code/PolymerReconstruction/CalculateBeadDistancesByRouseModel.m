classdef CalculateBeadDistancesByRouseModel<handle
    
    properties (Access=public)
        encounterMat
        encounterMatOriginal
        encounterProbOriginal
        connectivityMat
        params           % class for parameters
        graph
        chain           % the Rouse chain class
        smoother =Smoother; % signal smoother class
        nnEncounterProb = struct('distribution',[],'bins',[],'mean',[]);
        aboveLeft % probabilities above nn. encounter prob for left side
        aboveRight % probabilities above nn. encounter prob for right side 
    end
    
    properties (Access=private)

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
             
             obj.EstimateNearestNeighborEncounterProbability
             
             obj.AnalyzeConnectivity

            
            %  create a chain from the graph
            obj.CreateChainFromConnectivityGraph;
            obj.DisplayChain;
        end       
        
        function AnalyzeConnectivity(obj)
             % Perform analysis for each side 
             
             [eMatLeft, cMatLeft,obj.aboveLeft]    = obj.TransformEncountersToConnectivityMap('left');
             [eMatRight, cMatRight,obj.aboveRight] = obj.TransformEncountersToConnectivityMap('right');                                      
             eMat      = eMatLeft |eMatRight;
             obj.connectivityMat = eMat;
             cMat = cMatLeft |cMatRight;
             figure, imshow(cMat);
%             obj.DisplayConnectivityGraph(eMat,above,obj.distToAnalyze,obj.beadsToAnalyze);
        end
        
        function [eMat,cMat,above]= TransformEncountersToConnectivityMap(obj,side)
            % Analyze the encounters and determine bead neighbors in the
            % specified distance
            if strcmpi(side,'left')% left encounters
                % take the left side of the encounter matrix and rotate it 
                 encounters = obj.encounterMat(:,1:(size(obj.encounterMat,2)-1)/2);% tril(obj.encounterMat)';
                 encounters = fliplr(encounters);
            elseif strcmpi(side,'right')% right encounters 
                  encounters = obj.encounterMat(:,((size(obj.encounterMat,2)+1)/2 +1):end);
            end
            
            % Preallocations
            above = cell(1,size(encounters,1));% save indices of distances falling above the nearest neighor encounter probability
            dists = cell(size(encounters,1),size(encounters,2));
            histK = cell(size(encounters,1),size(encounters,2));
            chainRingStruct = cell(size(encounters,1),1);
            cMat  = false(size(obj.encounterMat,1)); 
            % Construct a binary connection matrix for a specific distance
            eMat = false(size(encounters,1),...
                         size(encounters,1),...
                         obj.params.reconstruction.numDistances);
            di          = diag(ones(1,size(eMat,1)-1),1)+diag(ones(1,size(eMat,1)-1),-1);% include nearest neighbors by default
            chainRingImage = zeros(size(encounters,1), size(encounters,2)*2 +1,3);  % an image representing positions of chains and rings                        
            connectivity = cell(size(encounters,1),1);
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
%                     [connectivity{bIdx}]          = obj.DetermineConnectivityByChainRingStruct(bIdx,chainRingStruct{bIdx},side);
                    
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
                    above{bIdx} = above{bIdx}(above{bIdx}>5); % don't allow the first 5 points to be included as peaks 
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
            
             
             for cIdx = 1:numel(connectivity)
                for c1Idx = 1:size(connectivity{cIdx})
                 cMat(connectivity{cIdx}(c1Idx,1), connectivity{cIdx}(c1Idx,2)) = true;
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
            [r,c]          = find(connectedBeads);    
            obj.params.chain.numBeads       = size(obj.encounterMat,1);
            obj.params.chain.connectedBeads = [r,c];   
             k = -(obj.params.chain.dimension* obj.params.chain.diffusionConst*...
                                                obj.params.chain.dt/obj.params.chain.b^2);
            obj.params.chain.springConst    = k*ones(size(obj.encounterMat,1)); 
            l  = logical(tril(obj.connectivityMat));
            u  = logical(triu(obj.connectivityMat));
            obj.connectivityMat = (l | l' | u | u');
            % set spring constants 
            for rIdx = 1:size(obj.connectivityMat,1)
                for cIdx = 1:size(obj.connectivityMat,1)
                    if (rIdx)>(cIdx) % left
                         above = obj.aboveLeft{(rIdx)};  
                        for aIdx = 1:numel(above)
                         if obj.connectivityMat((rIdx),(rIdx)-above(aIdx));
                             % ones(obj.params.chain.numBeads); % can be a scalar or a matrix the size of (numBeads) X (numBeads)
                          obj.params.chain.springConst((rIdx),(rIdx)-above(aIdx))=k*(obj.nnEncounterProb.mean/obj.encounterMat((rIdx),size(obj.encounterMat,1)-above(aIdx)))^1.5 *above(aIdx); 
                         end
                        end
                    elseif (cIdx)>(rIdx)% right
                        above = obj.aboveRight{(rIdx)};                        
                        for aIdx = 1:numel(above)
                           if obj.connectivityMat((rIdx),(rIdx)+above(aIdx));
                               % ones(obj.params.chain.numBeads); % can be a scalar or a matrix the size of (numBeads) X (numBeads)
                             obj.params.chain.springConst((rIdx),rIdx+above(aIdx))=k*(obj.nnEncounterProb.mean/obj.encounterMat((rIdx),size(obj.encounterMat,1)+1+above(aIdx)))^1.5 *above(aIdx); 
                           end
                        end
                    else % same bead 
                        obj.params.chain.springConst((rIdx),(cIdx))=-(obj.params.chain.dimension* obj.params.chain.diffusionConst*...
                                                obj.params.chain.dt/obj.params.chain.b^2);
                    end
                end
            end
            
            obj.chain = SimpleRouseFramework;
            obj.chain.Initialize(obj.params.chain);
%             obj.chain.Run;
            
        end
        
        function DisplayChain(obj)
            % Display the connected chain 
            if obj.params.chain.recordPath
                 ChainDynamicsPlayer(obj.chain.handles.classes.chain);
            else
                disp('encounters were not recorded, cannot display chain. Enable recordPath in the chain parameters')
            end
        end                    
        
        function ProcessEncounterMatrix(obj,encounterMat)
            % Interpolate and normalize the encounter histogram             
            obj.encounterMat                          = encounterMat;
            obj.encounterMatOriginal                  = encounterMat;
            obj.params.reconstruction.beadRange.bead1 = 1:size(encounterMat,1);
            obj.params.reconstruction.beadRange.bead2 = 1:size(encounterMat,2);
            
            obj.encounterMat(isnan(obj.encounterMat))= 0;% zero out nans

            % Get left and right parts of the encounter matrix 
            left  = fliplr(obj.encounterMat(:,1:(size(encounterMat,2)-1)/2));
            
            right = (obj.encounterMat(:,(size(encounterMat,2)+1)/2 +1 :end));
            
            % Save the original version of the prob after zero value
            % interpolation
            obj.encounterProbOriginal = [fliplr(left), zeros(size(left,1),1), right];
            
            % interpolate zero values
            [left,right] = obj.InterpolateZeroValues(left,right);                                  
            
            % smooth 'left' and 'right' parts 
            [left,right] = obj.CurveFitByInverseHeatEq(left,right);
            
            obj.encounterMat = [fliplr(left), zeros(size(left,1),1), right];
            obj.encounterMat(obj.encounterMat<0) = 0; % zero out negative values 
            
            % Normalize the smoothed encounter matrix
            for bIdx=obj.params.reconstruction.beadRange.bead1
                obj.encounterMat(bIdx,:) = InterpolateZeroValuesInSignal(obj.encounterMat(bIdx,:),obj.params.interpolation.zeroInterpolationMethod); % interpolate zero values
                obj.encounterMat(bIdx,size(obj.encounterMat,1))= 0;
                obj.encounterMat(bIdx,1:size(obj.encounterMat,1)-bIdx) = 0; % zero left size
                obj.encounterMat(bIdx,2*(size(obj.encounterMat,1))+1-bIdx:end) = 0;                
                obj.encounterMat(bIdx,:)= obj.encounterMat(bIdx,:)./obj.SumIgnoreNaN(obj.encounterMat(bIdx,:)); % normalize to get probabilities 
                obj.encounterMat(bIdx,isnan(obj.encounterMat(bIdx,:)))=0; % make nan equal zero
            end
        end
        
        function [left,right] = InterpolateZeroValues(obj,left,right)
            % Interpolate zero values in the signal                        
            for lIdx = 1:size(obj.encounterMat,1)
                left(lIdx,1:lIdx-1) = InterpolateZeroValuesInSignal(left(lIdx,1:lIdx-1),...
                                      obj.params.interpolation.zeroInterpolationMethod);
                left(lIdx,lIdx:end) = 0; % zero out the ends

            end
            
            right(isnan(right)) = 0;
            for rIdx = 1:size(right,1)
                right(rIdx,1:end-rIdx+1)     = InterpolateZeroValuesInSignal(right(rIdx,1:end-rIdx+1),...
                                            obj.params.interpolation.zeroInterpolationMethod);
                right(rIdx,end-rIdx+2:end)   = 0; 
            end 
            
            % Normalize the rows 
            s = sum(right+left,2);% normalization constant for each row
            for sIdx = 1:numel(s)
                if s(sIdx)>0
                    left(sIdx,:) = left(sIdx,:)./s(sIdx);
                    right(sIdx,:) = right(sIdx,:)./s(sIdx);
                end            
            end   
            
        end
        
        function [left,right] = CurveFitByInverseHeatEq(obj,left,right)
            % Smoothing parameters             
            regOrder    = 2;   % regularization order [0,1,2]
            lambda      = 2.7; % regularization constant
            alpha       = -1;
            numSpacePts = 3;
            initCond    = 3;
            minNumPts   = 5; % min number of points to perform analysis (smoothing)
                        
                        
            % perform smoothing
            for lIdx = minNumPts:size(obj.encounterMat,1)
               if ~all(left(lIdx,:)==0)% perform smoothing for non zero signals
                % smooth using the inverse heat equation                
                  [~,rL] = BoundaryElementHeatEquation('TestBemHeatEq_optimized',left(lIdx,1:lIdx-1),regOrder,lambda,alpha,numSpacePts,initCond,false); 
                  left(lIdx,1:lIdx-1) = rL./sum(rL);
               end
            end
            
            for rIdx =1:size(right)-minNumPts
               if ~all(right(rIdx,:)==0)
                  [~,rR] = BoundaryElementHeatEquation('TestBemHeatEq_optimized',right(rIdx,1:end-rIdx+1),regOrder,lambda,alpha,numSpacePts,initCond,false); 
                  right(rIdx,1:end-rIdx+1)=rR./sum(rR);
               end
            end
            
        end
        
        function EstimateNearestNeighborEncounterProbability(obj)
            % Use the encounter proability matrix to get an estimation of
            % the expected nearest neighor encounter probability 
            nnProb  = obj.encounterMat(:,[size(obj.encounterMat,1)-2:size(obj.encounterMat,1)-1 size(obj.encounterMat,1)+1:size(obj.encounterMat,1)+2]);
            nnProb = nnProb(nnProb>0);
            % calculate the histogram 
            [obj.nnEncounterProb.distribution,obj.nnEncounterProb.bins] = hist(nnProb(:),40);
             l = lognfit((nnProb+eps));
            obj.nnEncounterProb.mean = exp(l(1));
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
        
        function [connectivity] = DetermineConnectivityByChainRingStruct(obj,beadNumber, chainRingStruct,side)
            % Determine which beads are connected according to the results
            % in the chainRingStruct. Connected beads are those at the end
            % of rings
            connectivity = [];
            
            for cIdx = 1:numel(chainRingStruct)
                if strcmpi(chainRingStruct(cIdx).type,'ring')
                    if strcmpi(side,'left')
                        connectivity(end+1,:) = beadNumber-[chainRingStruct(cIdx).endInd chainRingStruct(cIdx).startInd];
                    elseif strcmpi(side,'right')
                        connectivity(end+1,:) = beadNumber+[chainRingStruct(cIdx).startInd chainRingStruct(cIdx).endInd];
                    else
                        error('unsupported side option')
                    end
                end           
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
%                 modelValues = modelValues./modelValues(1) *max(prob(1:min(15,size(prob,2))));
                dists       = (prob *s).^(-1./beta);
            elseif strcmpi(obj.params.reconstruction.prob2distMethod,'rouse')
                % use rouse with beta=1.5
                beta        = 1.5;
                inds        = find(~isnan(prob));
                s           = sum(inds.^(-beta));
                modelValues = obj.params.optimization.model(beta,inds);
                modelValues = modelValues./modelValues(1) * obj.nnEncounterProb.mean;%max(prob(1:min(15,size(prob,2))));
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
        
        function [observedProb] = ProcessBeadEncounterSignal(obj,encounterSignal)%obsolete
            % process a single encounter signal
            % fill NaN position with nearest neighbors mean value
            observedProb = InterpolateZeroValuesInSignal(encounterSignal,obj.params.interpolation.zeroInterpolationMethod);
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