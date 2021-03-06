classdef SimpleRouseFrameworkNew<handle
    properties
        params
        beadDistance
        encounterHistogram
        beadEncounterHistogram
        beadEncounterProbability
        meanEncounterProbability
        fitModel 
        goalFunction
        fitResults
        fOptions% fit options
        handles
        recipe % actions defined in the recipe file 
        round
        simulation = 0;
        step       = 0;
        simulationData
        stepExitFlag = false;
        userData% various user data needed for passing between functions
    end
    
    events
    end
    
    methods
        function obj = SimpleRouseFramework
        end
        
        function Initialize(obj,rouseParams)
            % set input parameters
            obj.params              = rouseParams;
            
            % set fit options 
            obj.fitModel             = @(b,dists)(1/sum(dists.^(-b)))*dists.^(-b);
            obj.goalFunction         = @(b)sum((obj.fitModel(b,1:obj.params.numBeads-1)-obj.fitModel(1.5,1:obj.params.numBeads-1)).^2);
            obj.fOptions             = optimset; 
            obj.fOptions.TolX        = 1e-16;
            obj.fOptions.TolFun      = 1e-16;
            obj.fOptions.MaxFunEvals = 2500;
            obj.fOptions.MaxIter     = 2500;
                       
            obj.round = 0; % start with zero index, the values increases in PreRoundActions function
            % preallocations
            obj.fitResults.mean           = cell(1,obj.params.numRounds);
            obj.fitResults.bead.fittedExp = zeros(obj.params.numBeads,obj.params.numRounds);
            obj.fitResults.bead.gof       = cell(obj.params.numBeads,obj.params.numRounds);
            obj.beadEncounterHistogram    = zeros(obj.params.numBeads,obj.params.numBeads-1,obj.params.numRounds);
            obj.encounterHistogram        = zeros(obj.params.numBeads,obj.params.numBeads,obj.params.numRounds);  
            obj.beadEncounterHistogram    = zeros(obj.params.numBeads,obj.params.numBeads-1,obj.params.numRounds);
            obj.beadEncounterProbability  = zeros(obj.params.numBeads,obj.params.numBeads-1,obj.params.numRounds);
            obj.meanEncounterProbability  = zeros(1,obj.params.numBeads-1,obj.params.numRounds); 
            
             % Initialize sequence 
            obj.ReadRecipeFile;
            obj.Run;
            obj.AnalyzeResults;   
            obj.DisplayAnalyzedData
        end
        
        function ReadRecipeFile(obj)
            
         t = fileread(fullfile(obj.params.recipeFolder, [obj.params.recipeFileName,'.rcp']));
         % search for the function marker 
         [funcStartPos1,funcStartPos2] = regexp(t,'<func>');
         [funcEndPos1,funcEndPos2]     = regexp(t,'</func>');
         for fIdx = 1:numel(funcStartPos1)-1
             % sort the functions into categories
             functionName = strtrim(t(funcStartPos2(fIdx)+1:funcEndPos1(fIdx)-1));
             obj.recipe.(functionName) = t(funcEndPos2(fIdx)+1:funcStartPos1(fIdx+1)-1);
         end
             functionName = strtrim(t(funcStartPos2(end)+1:funcEndPos1(end)-1));
             obj.recipe.(functionName) = t(funcEndPos2(end)+1:end);     
        end
        
        function PreRoundActions(obj)
             obj.round = obj.round+1;
             eval(obj.recipe.PreSimulationBatchActions);
             
        end
        
        function PreRunActions(obj)
            % actions performed before each simulation run
            obj.simulation = obj.simulation+1;
            eval(obj.recipe.PreRunActions);
        end
        
        function PreStepActions(obj)
            % actions before each simulation step
          obj.step = obj.step+1;
          eval(obj.recipe.PreStepActions);
        end
        
        function Run(obj)
            % run simulations
            % preallocation 
            obj.beadDistance = zeros(obj.params.numBeads,obj.params.numBeads,obj.params.numSimulations,obj.params.numRounds);
            
            for rIdx = 1:obj.params.numRounds
                obj.PreRoundActions;
                tic
              for sIdx = 1:obj.params.numSimulations
                obj.PreRunActions
                obj.handles.classes.rouseChain = SimpleRouse(obj.params);
                obj.handles.classes.rouseChain.Initialize
                while obj.step<=obj.params.numSteps && ~obj.stepExitFlag
                 obj.PreStepActions
                 obj.handles.classes.rouseChain.Step;
                 obj.PostStepActions
                end
%                 obj.handles.classes.rouseChain.Run;% run rouse
                obj.beadDistance(:,:,sIdx,rIdx) = obj.handles.classes.rouseChain.beadDist(:,:,end);% take the last recorded bead distances
               
                obj.simulationData(obj.round).simulation(sIdx).time = obj.handles.classes.rouseChain.simulationTime;
                obj.simulationData(obj.round).simulation(sIdx).numSteps = obj.handles.classes.rouseChain.params.numSteps;
                obj.simulationData(obj.round).simulation(sIdx).meanStepTime = obj.handles.classes.rouseChain.simulationTime/obj.handles.classes.rouseChain.params.numSteps;
                
                obj.PostRunActions;% perform post run actions
                sprintf('%s%d%s%d%s', 'round ', rIdx,' simulation ', sIdx, ' is done')
              end
              
              % record simulation data
              obj.simulationData(obj.round).round.time = toc;% record round time 
              obj.simulationData(obj.round).params     = obj.params;% save round parameters
              
              % perform user defined actions at the end of each round 
              obj.PostRoundActions
            end
        end
        
        function PostStepActions(obj)
         % actions performed after each step 
         eval(obj.recipe.PostStepActions);
        end
        
        function PostRunActions(obj)
          % actions performed after each simulation
          eval(obj.recipe.PostRunActions)
          
          obj.step =0; % restart step index;
        end
        
        function PostRoundActions(obj)
            % perform an action after the simulation round is over 
             eval(obj.recipe.PostSimulationBatchActions); 
             obj.simulation = 0; % reset the simulation counter
        end
        
        function AnalyzeResults(obj)
            obj.CalculateEncounterHistogram;
            obj.CalculateBeadEncounterHistogram;
            obj.CalculateBeadEncounterProbability;
            obj.CalculateMeanEncounterProbability;
            obj.FitMeanEncounterProbability
            obj.FitBeadEncounterProbability;

        end
        
        function DisplayAnalyzedData(obj)
            obj.DisplayFittedExponents;
            obj.DisplayMeanEncounterProbability
            obj.DisplayBeadDataFit
            obj.DisplayEncounterHistograms
        end
        
        function CalculateEncounterHistogram(obj)
            % Calculate the encounter histogram over all rounds 
            for rIdx = 1:obj.params.numRounds
             % Sum over all experiments
             obj.encounterHistogram(:,:,rIdx) = sum(obj.beadDistance(:,:,:,rIdx)<obj.params.encounterDist,3);%  encounters 
             % remove main diagonal 
             obj.encounterHistogram(:,:,rIdx) = obj.encounterHistogram(:,:,rIdx)-diag(diag(obj.encounterHistogram(:,:,rIdx)));
            end
        end
        
        function CalculateBeadEncounterHistogram(obj)
            % Calculate the encounter histogram for each bead in each round          
            for rIdx = 1:obj.params.numRounds                
                for bIdx =1:obj.params.numBeads
                    k1 = numel(1:numel(bIdx+1:obj.params.numBeads));
                    k2 = numel(1:numel(1:bIdx-1));             
                    
                    f = zeros(2,obj.params.numBeads-1); % sum 'right' and 'left' encounters
                    f(1,1:numel(bIdx+1:obj.params.numBeads)) = obj.encounterHistogram(bIdx,bIdx+1:obj.params.numBeads,rIdx);% right side
                    f(2,1:numel(1:bIdx-1)) = fliplr(obj.encounterHistogram(bIdx,1:bIdx-1,rIdx));% flipped left side
                    f    = sum(f);% sum left and right

                    normFac = ones(1,size(f,2));
                    normFac(1:min([k1,k2]))=2;
                    f = f./normFac;
                    obj.beadEncounterHistogram(bIdx,:,rIdx)=f;
                end
                
            end
        end
        
        function CalculateBeadEncounterProbability(obj)
            % calculate the bead encounter probability by dividing the
            % beadEncounterHistogram by the sumof its rows            
%             obj.beadEncounterProbability = obj.beadEncounterHistogram;    
            % divide by the sum of each row to get the probability
            for pIdx = 1:obj.params.numRounds
                for bIdx = 1:obj.params.numBeads
                    if sum(obj.beadEncounterHistogram(bIdx,:,pIdx))~=0
                        % divide by the sum of the histogram row to get
                        % probability 
                        obj.beadEncounterProbability(bIdx,:,pIdx) =obj.beadEncounterHistogram(bIdx,:,pIdx)./sum(obj.beadEncounterHistogram(bIdx,:,pIdx));
                    else
                        obj.beadEncounterProbability(bIdx,:,pIdx) =obj.beadEncounterHistogram(bIdx,:,pIdx);
                    end
                end
            end
        end
        
        function CalculateMeanEncounterProbability(obj)
                % Calculate mean encounter probability
                % the normalization term is determined by the number of terms avialable
                % for each bead.
            for rIdx = 1:obj.params.numRounds
                normTerm = obj.params.numBeads*ones(1,obj.params.numBeads-1);
                normTerm(1:floor(obj.params.numBeads/2)) = obj.params.numBeads;
                s = obj.params.numBeads-2:-1:1;
                % take any second element of s
                s= s(1:2:end);                
                normTerm(1,end-numel(s)+1:end)=s; 
                obj.meanEncounterProbability(1,:,rIdx) = sum(obj.beadEncounterHistogram(:,:,rIdx),1)./normTerm;
                obj.meanEncounterProbability(1,:,rIdx) = obj.meanEncounterProbability(1,:,rIdx)./sum(obj.meanEncounterProbability(1,:,rIdx));
            end
        end        
        
        function FitMeanEncounterProbability(obj)
            % fit the theoretical mean encounter curve for all beads int he
            % current round. Reminder, the round includes several
            % simulations 

            dists  = (1:obj.params.numBeads-1)';
            for rIdx = 1:obj.params.numRounds
%                 obj.fOptions.Weights = obj.meanEncounterProbability(1,:,pIdx)~=0;
%                 [obj.fitResults.mean{pIdx}, ~] = fit(dists,obj.meanEncounterProbability(1,:,pIdx)',obj.fitModel,...
%                     obj.fOptions);
                f =@(b)sum(((1/sum((dists).^(-b)))*(dists).^(-b)-obj.meanEncounterProbability(1,:,rIdx)').^2);
                [obj.fitResults.mean{rIdx}.b] = fminbnd(f,0.2,2.2,obj.fOptions);
%                 obj.fOptions.Weights = [];
            end
          
            
        end        
        
        function FitBeadEncounterProbability(obj)
            % fit the encounter probability curve for each bead in each
            % round.Reminder: a round includes several simulations 
            
            dists = (1:obj.params.numBeads-1)';
            for pIdx = 1:obj.params.numRounds
             for bIdx = 1:obj.params.numBeads
                 k1 = numel(1:numel(bIdx+1:obj.params.numBeads));
                 k2 = numel(1:numel(1:bIdx-1));    
                 minK = min([k1,k2]);
                 numPts = obj.params.numBeads-1-minK;% take as many points as there is data available (omit lagging zeros)
%                 [beadFit, gof] = fit(dists,obj.beadEncounterProbability(bIdx,:,pIdx)',obj.fitModel,obj.fOptions);
                f = @(b)sum((((1/sum((dists(1:numPts)).^(-b)))*(dists(1:numPts)).^(-b))-obj.beadEncounterProbability(bIdx,1:numPts,pIdx)').^2);
                e = fminbnd(f,0.2,2.2,obj.fOptions);

                obj.fitResults.bead.fittedExp(bIdx,pIdx) = e;%beadFit.b;
%                 obj.fitResults.bead.fittedBias(bIdx,pIdx) = beadFit.a;
%                 obj.fitResults.bead.gof{bIdx,pIdx} = gof;
             end
            end
        end
        
        function DisplayFittedExponents(obj)
            figure, 
            plot(obj.fitResults.bead.fittedExp,'.');
            xlabel('Bead Number','FontSize',16);
            ylabel('Fitted exponent','Fontsize',16);
        end
        
        function DisplayBeadDataFit(obj,beadRange)
            if ~exist('beadRange','var')
                beadRange = 1:obj.params.numBeads;                
            end
            
            for rIdx = 1:obj.params.numRounds
            f = figure('FileName',sprintf('%s%s','BeadDataFit_Experiment',num2str(rIdx)));
            a = axes('Parent',f,'NextPlot','Add');
            dists = 1:obj.params.numBeads-1;
            for bIdx = 1:numel(beadRange)
                lineColor = rand(1,3);
                % plot the Encounter probability from simulation 
                line('Parent',a,...
                     'XData',dists,...
                     'YData',obj.beadEncounterProbability(bIdx,:,rIdx),...
                     'Color',lineColor,...
                     'DisplayName',sprintf('%s%s','Bead',num2str(bIdx)));
                 % plot the fitted line 
                 line('Parent',a,...
                     'XData',dists,...
                     'YData',(1/sum(dists.^(-obj.fitResults.bead.fittedExp(bIdx,rIdx))))*dists.^(-obj.fitResults.bead.fittedExp(bIdx,rIdx)),...%obj.fitModel(obj.fitResults.bead.fittedBias(bIdx,rIdx),obj.fitResults.bead.fittedExp(bIdx,rIdx),dists),...
                     'Color',lineColor,...
                     'DisplayName',sprintf('%s%s','exp= ',num2str(obj.fitResults.bead.fittedExp(bIdx,rIdx))));                                  
             end
            end
            
            % display expected model 
            line('XData',dists,...
                'YData',(1/sum(dists.^(-1.5)))*dists.^(-1.5),...%obj.fitModel(-0.5/(-1+obj.params.numBeads^-0.5),1.5,dists),...
                'LineWidth',3,...
                'Color','r',...
                'Parent',a,...
                'DisplayName','Theoretical model')
            
            xlabel('Bead distance','Fontsize',16);
            ylabel('Encounter prob','FontSize',16);
        end
        
        function DisplayMeanEncounterProbability(obj)
            f = figure('FileName','meanEncounterProbability') ;
            a = axes('Parent',f,'NextPlot','Add');
            dists = (1:obj.params.numBeads-1)';
            for rIdx = 1:obj.params.numRounds
                c = rand(1,3);
                line('XData',dists,...
                    'YData',obj.meanEncounterProbability(1,:,rIdx),...
                    'Color',c,...                  
                    'Parent',a),
                
                line('XData',dists,...
                    'YData',(1/sum(dists.^(-1.5)))*dists.^(-1.5),...%obj.fitModel(obj.fitResults.mean{rIdx}.a,obj.fitResults.mean{rIdx}.b,dists),...
                    'Color',c,...
                    'DisplayName',[ 'b=', num2str(obj.fitResults.mean{rIdx}.b)],...
                    'Parent',a)
%                 m(rIdx) = obj.fitResults.mean{rIdx}.a;
            end
            % add the theoretical model 
            line('XData',dists,...
                'YData',(1/sum(dists.^(-1.5)))*dists.^(-1.5),...
                'Color','r',...
                'Parent',a,...
                'DisplayName','theoretical encounter curve')
        end
        
        function DisplayEncounterHistograms(obj)
            
            for rIdx = 1:obj.params.numRounds
                figure, imagesc(obj.encounterHistogram(:,:,rIdx)), colormap hot
                title(sprintf('%s%d','experiment ', rIdx),'FontSize',16);
                xlabel('Bead number','FontSize',16);
                ylabel('BeadNumber','FontSize', 16);
            end
            
        end
    end
end