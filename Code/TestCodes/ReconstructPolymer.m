classdef ReconstructPolymer<handle
    
    properties
        encounterData
        encounterProb
        firstNeighborProbMean
        firstNeighborProbStd
        firstNeighborProbDist
        firstNeighborWbl
    end
    
    methods        
        
        function obj = ReconstructPolymer
        end
        
        function Initialize(obj)
            load('savedAnalysisTADDAndE');
            % this class relies on the presence of the saved analysis 
            obj.encounterData = a.beadData.encounterData.twoSides.average;
            obj.EstimateClosestNeighborEncounterProbability;
        end
        
        function EstimateClosestNeighborEncounterProbability(obj)
            f              = zeros(size(obj.encounterData,1),2);
            fEncounterProb = zeros(size(obj.encounterData,1),1);
            for eIdx = 1:size(obj.encounterData,1)
                obj.encounterProb(eIdx,:) = [fliplr(obj.encounterData{eIdx,1}) obj.encounterData{eIdx,2}];
                sumSig = sum(obj.encounterProb(eIdx,:));
                 obj.encounterProb(eIdx,:) = obj.encounterProb(eIdx,:)/sumSig;
                if ~isempty(obj.encounterData{eIdx,1})
                  f(eIdx,1) = obj.encounterData{eIdx,1}(1)./sumSig;
                end                
                if ~isempty(obj.encounterData{eIdx,2})
                  f(eIdx,2) = obj.encounterData{eIdx,2}(1)./sumSig;
                end
                
                if ~any(isnan(f(eIdx,:)))
                    zInds = ~(f(eIdx,:)==0);
                    if any(zInds)
                      fEncounterProb(eIdx) = mean(f(eIdx,zInds));
                    end
                end
            end
            % fit Weibull distribution to first neighbor encounter Prob and
            % extract mean and std
            obj.FitWeibullToEncounterDataOfFirstNeighbor(fEncounterProb);                        
        end
        
        function FitWeibullToEncounterDataOfFirstNeighbor(obj,fEncounterProb)
            % Fit a Weibull distribution to the encounter data of the first
            % neighbor 
            cens = fEncounterProb==0;
            fEncounterProb(cens) = eps;
            [obj.firstNeighborWbl] = fitdist(fEncounterProb,'wbl','Censoring',cens);
            lambda = obj.firstNeighborWbl.A;
            k      = obj.firstNeighborWbl.B;
            obj.firstNeighborProbMean = lambda*gamma(1+1/k);
            obj.firstNeighborProbStd  = sqrt((lambda^2)*(gamma(1+2/k)-gamma(1+1/k)^2));
            obj.firstNeighborProbDist = fEncounterProb;
        end
        
        function DisplayEncounterProb(obj,beadNum)
            f = figure('Name',sprintf('%s%s','encounter Prob'));
            a = axes('Parent',f,'Linewidth',3,'FontSize',30);
            for bIdx = 1:numel(beadNum)
                lineC = rand(1,3);
                line('XData',[-beadNum(bIdx)+1:1:-1, 1:307-beadNum(bIdx)],...
                    'YData',obj.encounterProb(beadNum(bIdx),:),...
                    'Parent',a,...
                    'Color',lineC,...
                    'lineWidth',4,...
                    'DisplayName',(sprintf('%s%s','encounter Prob bead ',num2str(beadNum(bIdx)))))
            end
            axis tight
            % plot the mean encounter prob and std   
            xLim = get(a,'XLim');
            line('XData',xLim,...
                 'YData',[obj.firstNeighborProbMean, obj.firstNeighborProbMean],...
                 'Parent',a,...
                 'DisplayName','\mu');
            line('XData',xLim,...
                 'YData',[obj.firstNeighborProbMean-obj.firstNeighborProbStd, obj.firstNeighborProbMean-obj.firstNeighborProbStd],...
                 'Parent',a,...
                 'LineStyle','-.',...
                 'DisplayName','\mu-\sigma');
            line('XData',xLim,...
                 'YData',[obj.firstNeighborProbMean+obj.firstNeighborProbStd, obj.firstNeighborProbMean+obj.firstNeighborProbStd],...
                 'Parent',a,...
                 'LineStyle','-.',...
                 'DisplayName','\mu+\sigma');
             
            axis tight
%             legend(get(a,'Children'));
        end
        
    end
end