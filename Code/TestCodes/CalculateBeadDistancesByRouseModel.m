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
        graph
        model
        fitOpt
        dataFolder
        dataFileName
    end
    
    methods
        
        function obj = CalculateBeadDistancesByRouseModel
        end
        
        function SetDefaultParams(obj)
            obj.beadRange      = struct('bead1',1:70,...
                                        'bead2',1:70);
            obj.smoothingSpan  = 1;
            obj.smoothingMethod= 'loess'; % see smooth function for options 
            obj.numDistances   = 20;       % for how many distances to perform analysis for connectivity
            obj.distToAnalyze  = 4;       % can be a vector of integers, for what disance to show the analysis
            obj.beadsToAnalyze = 5;   % for what beads to show the connectivity graphs
            obj.model          = fittype('(1/sum(x.^(-beta))).*x.^(-beta)');
            obj.fitOpt         = fitoptions(obj.model);
            set(obj.fitOpt,'Lower',0,'Upper',1.5,'StartPoint',1,'Robust','off');
            obj.dataFolder     = 'D:\Ofir\Work\ENS\Polymer Chain Dynamics\Code\ExperimentDataAnalysis';
            obj.dataFileName   = 'savedAnalysisTADDAndE';
        end
        
        function Initialize(obj)
            obj.SetDefaultParams;
            load(fullfile(obj.dataFolder,obj.dataFileName))
            [~,~,obj.encounterMat,~] = a.ProcessEncounters(obj.beadRange,'average');
            % Truncate the encounter matrix
            obj.encounterMat = obj.encounterMat(obj.beadRange.bead1,obj.beadRange.bead2(1:(end-1))-obj.beadRange.bead1(1)+1);
            
            % preallocations
            above = cell(1,numel(obj.beadRange.bead1));% save indices of distances falling above the nearest neighor encounter probability
            dists = cell(size(obj.encounterMat,1),size(obj.encounterMat,2));
            histK = cell(size(obj.encounterMat,1),size(obj.encounterMat,2));
            beta  = zeros(size(obj.encounterMat,1),1);
            % construct a binary connection matrix for a specific distance
            eMat = false(numel(obj.beadRange.bead1),numel(obj.beadRange.bead2),obj.numDistances);
            di   = eMat(:,:,1);%diag(ones(1,size(eMat,2)-1),1)+diag(ones(1,size(eMat,2)-1),-1);% include nearest neighbors by default
            
            
            for bIdx = 1:size(obj.encounterMat,1);
                
                observedProb = obj.SmoothSignal(obj.encounterMat(bIdx,~isnan(obj.encounterMat(bIdx,:))),obj.smoothingSpan,obj.smoothingMethod)';
                
                if ~all(isnan(observedProb))
                    
                    observedProb = observedProb/sum(observedProb);
                    
                    % Divide the probabilites into distances according to a division given by the
                    % expected model
                    inds        = find(~isnan(observedProb));
                     [fitStruct] = fit(inds',observedProb(inds)',obj.model,obj.fitOpt);
                    beta(bIdx)  = fitStruct.beta;
                    k           =  obj.model(beta(bIdx),inds);
                    % normalize to match the nearest neighbor encounter probability
                    if mod(bIdx,5)==0
                        obj.PlotBeadClusteringByDistance(observedProb,inds,k);
                        title(num2str(bIdx))
                    end
                    % Calculate the histogram
                    above{bIdx} = find(observedProb>k(1));
                    below{bIdx} = [];
                    for kIdx = 1:numel(k)-1
                        dists{bIdx,kIdx} = find(observedProb>k(kIdx+1) & observedProb<=k(kIdx));
                        below{bIdx,kIdx} = dists{bIdx,kIdx}((dists{bIdx,kIdx}<kIdx));                                                
                        % add the terms "above" to the dist 1 neighbors
                        histK{bIdx,kIdx} = numel(dists{bIdx,kIdx});
                        if kIdx ==1
                            dists{bIdx,kIdx} = [dists{bIdx,kIdx} above{bIdx}];
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
            obj.DisplayConnectivityGraph(eMat,above,obj.distToAnalyze,obj.beadsToAnalyze);
            
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
            inds                        = setdiff(1:size(eMat,1),beadToAnalyze);
            obj.connectivityMat(inds,:) = false;
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
    
    methods (Static)
        
        function m = MeanIgnoreNaN(sigIn)
            sigIn(isnan(sigIn)) = 0;
            m = mean(sigIn);
            m = m./sum(m);
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
        
        function sigOut = SmoothSignal(sigIn,smoothSpan,method)
            % Smooth a signal sigIn with a smoothing span smoothSpan
            sigOut = smooth(sigIn,smoothSpan,method);
        end
    end
    
end