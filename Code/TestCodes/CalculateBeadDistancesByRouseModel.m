classdef CalculateBeadDistancesByRouseModel<handle
    
    properties
        encounterMat
        connectivityMat
        beadRange       % set bead range for TAD D
        smoothingMethod % see smooth function for options 
        smoothingSpan   % smoothing span for the encounter probability signal
        numDistances    % for how many distances to perform analysis for connectivity
        distToAnalyze   % can be a vector of integers, for what disance to show the analysis
        beadsToAnalyze  % for what beads to show the connectivity graphs
        zeroInterpolationMethod % how to fill in nan values of the encounter signal [ linear, nearest, spline, cubic,...]
        graph
        model
        chain           % the Rouse chain class
        fitOpt
        dataFolder
        dataFileName
    end
    properties (Access=private)
        smoother = Smoother;
    end
    
    methods
        
        function obj = CalculateBeadDistancesByRouseModel
        end
        
        function SetDefaultParams(obj)
            obj.beadRange      = struct('bead1',1:307,...
                                        'bead2',1:307);
            obj.smoothingSpan  = 15;      % can be a vector of integer, in the present version only the last span is considered
            obj.smoothingMethod= 'gaussian'; % see smooth function for options 
            obj.numDistances   = 1;       % for how many distances to perform analysis for connectivity
            obj.distToAnalyze  = [1];     % can be a vector of integers, for what disance to show the analysis
            obj.beadsToAnalyze = 1;      % for what beads to show the connectivity graphs
            obj.zeroInterpolationMethod = 'linear'; % how to fill in gaps for zero values of the encoutner signal (see interp1 documentation for methods)
            obj.model          = fittype('(1/sum(x.^(-beta))).*x.^(-beta)');
            obj.fitOpt         = fitoptions(obj.model);
            set(obj.fitOpt,'Lower',0,'Upper',1.5,'StartPoint',1,'Robust','off');
            obj.dataFolder     = fullfile(pwd,'ExperimentDataAnalysis');
            obj.dataFileName   = 'savedAnalysisTADDAndE';
        end
        
        function Initialize(obj,encounterMat)
            obj.SetDefaultParams;
            if ~exist('encounterMat','var')
            load(fullfile(obj.dataFolder,obj.dataFileName))
            [~,~,obj.encounterMat,~] = a.ProcessEncounters(obj.beadRange,'average');
            
            else
                obj.encounterMat    = encounterMat;
                obj.beadRange.bead1 = 1:size(encounterMat,1);
                obj.beadRange.bead2 = 1:size(encounterMat,2);
            end
            % Truncate the encounter matrix according to the bead range specified 
            obj.encounterMat = obj.encounterMat(obj.beadRange.bead1,obj.beadRange.bead2-obj.beadRange.bead1(1)+1);
            
            % preallocations
            above = cell(1,numel(obj.beadRange.bead1));% save indices of distances falling above the nearest neighor encounter probability
            dists = cell(size(obj.encounterMat,1),size(obj.encounterMat,2));
            histK = cell(size(obj.encounterMat,1),size(obj.encounterMat,2));
            beta  = zeros(size(obj.encounterMat,1),1);
            % construct a binary connection matrix for a specific distance
            eMat = false(numel(obj.beadRange.bead1),numel(obj.beadRange.bead1),obj.numDistances);
            di   = diag(ones(1,size(eMat,1)-1),1)+diag(ones(1,size(eMat,1)-1),-1);% include nearest neighbors by default
            
            % get expected signal for all encounters 
            [expectedSignal,fitStructModel] = obj.GetExpectedEncounterSignal(obj.encounterMat);
                        
            for bIdx = 1:size(obj.encounterMat,1);% for each bead
                
                observedProb = obj.ProcessBeadEncounterSignal(obj.encounterMat(bIdx,:));                
                                
                if ~all(isnan(observedProb))                                        
%                     expectedSignal = expectedSignal*sum(expectedSignal)/sum(observedProb);
                    % Divide the probabilites into distances according to a division given by the
                    % expected model
                    inds        = find(~isnan(observedProb));
                    [fitStruct] = fit(inds',observedProb(inds)',obj.model,obj.fitOpt);
                    beta(bIdx)  = fitStruct.beta;
                    s           = sum(inds.^(-beta(bIdx)));
                    modelValues = obj.model(beta(bIdx),inds);
                    modelValues = modelValues./modelValues(1) *max(observedProb(1:10));
%                     observedProb= observedProb*modelValues(1)/observedProb(1);

                    % normalize to match the nearest neighbor encounter probability
                    if mod(bIdx,408)==0
                        obj.PlotBeadClusteringByDistance(observedProb,inds,modelValues);
                        title(num2str(bIdx))
                    end
                    % Calculate the histogram
                    above{bIdx} = find(observedProb>modelValues(1));
                    below{bIdx} = [];
                    for kIdx = 1:numel(modelValues)-1
                        dists{bIdx,kIdx} = find(observedProb>modelValues(kIdx+1) & observedProb<=modelValues(kIdx));
                        below{bIdx,kIdx} = dists{bIdx,kIdx}((dists{bIdx,kIdx}<kIdx));                                                
                        % add the terms "above" to the dist 1 neighbors
                        histK{bIdx,kIdx} = numel(dists{bIdx,kIdx});
                        if kIdx ==1
                            dists{bIdx,kIdx} = [dists{bIdx,kIdx}, above{bIdx}];
                        end
                    end
                    
                end
            end
            
            for dIdx = obj.distToAnalyze
                eMat(:,:,dIdx) = eMat(:,:,dIdx)|di;
                for b1Idx = 1:size(obj.encounterMat,1)
                    % collect all beads at distance 1
                    inds1 = b1Idx+dists{b1Idx,dIdx};
                    inds1 = inds1(inds1<numel(obj.beadRange.bead2));
                    inds2 = b1Idx-dists{b1Idx,dIdx};
                    inds2 = inds2(inds2>=1);
                    eMat(b1Idx,[inds1 inds2],dIdx)= true;
                end
            end
            
            obj.connectivityMat = eMat;
%             obj.DisplayConnectivityGraph(eMat,above,obj.distToAnalyze,obj.beadsToAnalyze);
            
            %  create a chain from the graph
            obj.CreateChainFromConnectivityGraph;
            obj.DisplayChain;
        end
        
        function [expectedSignal,fitStructModel] = GetExpectedEncounterSignal(obj,encounterMat)
            % Get the expected encounter signal from an encounter matrix 
            expectedSignal    = obj.MeanIgnoreNaN(encounterMat);
            expectedSignal    = obj.smoother.Smooth(expectedSignal',obj.smoothingMethod,obj.smoothingSpan);
            se                = obj.SumIgnoreNaN(expectedSignal(:,1,1));% normalize
            expectedSignal    = expectedSignal(:,1,1)./se;
            [fitStructModel]  = fit((1:numel(expectedSignal))',expectedSignal,obj.model,obj.fitOpt);
        end
        
        function [observedProb] = ProcessBeadEncounterSignal(obj,encounterSignal)
            % process a single encounter signal 
            % fill NaN position with nearest neighbors mean value
            observedProb = obj.InterpolateZeroValuesInSignal(encounterSignal);
            observedProb = obj.smoother.Smooth(observedProb,obj.smoothingMethod,obj.smoothingSpan);
%             observedProb = obj.SmoothSignal(observedProb,...
%                                             obj.smoothingSpan,obj.smoothingMethod);
            sop            = obj.SumIgnoreNaN(observedProb(1,:,end));
            observedProb   = observedProb(1,:,end)./sop; % normalize
        end
        
        function GetDistanceDistribution(obj,dist)
            % calculate the distribution for a specifind distance
            figure,hold on
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
                        'LineWidth',4);
                end
                
            end
            legend(get(gca,'Children'));
            
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
            connectedBeads = obj.connectivityMat;
            % remove the trivial connections on the super diagonal        
            connectedBeads= triu(connectedBeads-diag(diag(connectedBeads,1),1));          
            [r,c]       = find(connectedBeads);    
            sr          = SimpleRouseParams;
            sr.numBeads = size(obj.encounterMat,1);                       
            sr.recordPath     = true;
            sr.numBeads       = numel(obj.beadRange.bead1);
            sr.dt             = 1e-2;
            sr.numSteps       = 500;
            sr.noiseSTD        = 0.00;
            sr.dimension       = 3;      
            sr.b               = sqrt(3);
            sr.diffusionConst  = 1;
            sr.numSimulations  = 1;
            sr.springConst     = -(sr.dimension* sr.diffusionConst* sr.dt/ sr.b^2)*ones( sr.numBeads); % can be a scalar or a matrix the size of (numBeads) X (numBeads)
            sr.connectedBeads  = [r c];
            obj.chain = SimpleRouse(sr);
            obj.chain.Initialize;
            obj.chain.Run;
            
        end
        
        function DisplayChain(obj)
            % display the connected chain 
            ChainDynamicsPlayer(obj.chain);
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
        
    end
    
    methods (Access=private)
        function [sigOut] = InterpolateZeroValuesInSignal(obj,sigIn)
            % remove nan values by interpolation of the signal 
            sigOut     = sigIn;
            zeroInds   = find(sigIn==0);
            noZeroInds = find(~(sigIn==0) & ~isnan(sigIn));
            if ~isempty(zeroInds)
                % interpolate the signal in the nan positions 
            x     = noZeroInds;
            y     = sigIn(x);
            sigOut(zeroInds)= interp1(x,y,zeroInds,obj.zeroInterpolationMethod);
            end            
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
        
        function sigOut = SmoothSignal(sigIn,smoothSpan,method)%obsolete
            % Smooth a signal sigIn with a smoothing span smoothSpan
            % try all smoothing spans defing in the parameter smoothiSpan
            ind    = 1;
            sigOut = zeros(1,numel(sigIn),numel(smoothSpan));
            for s = smoothSpan
            sigOut(1,:,ind) = smooth(sigIn,s,method);
            ind = ind+1;
            end
        end
    end
    
end