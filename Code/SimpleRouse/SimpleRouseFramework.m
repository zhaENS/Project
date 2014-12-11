classdef SimpleRouseFramework<handle
    % this is a framework to run and analyze the results of a simple rouse
    % chain 
    %TODO: finish the connection by Brownian bridge function 
    
    properties
        params
        beadDistance
        msd % values of <|x(t)-x(0)|^2> with x(t) the position of each bead
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
        polymerModel %string 'betaChain' or 'RouseChain'
        userData% various user data needed for passing between functions
    end
    
     properties (Access = private)   
     end     
     
    events
    end
    
    methods
        function obj = SimpleRouseFramework
        end
        
        function Initialize(obj,rouseParams)
            % set input parameters
            obj.params  = rouseParams;
            obj.round   = 0; % start with zero index, the values increases in PreRoundActions function

            
            % Initialize sequence
            obj.ReadRecipeFile;   % read recipe file 
            obj.SetRecipeParams;  % load the parameters defined in the recipe file 
            if obj.params.beta~=2
                obj.polymerModel = 'betaChain';
            else 
                obj.polymerModel = 'rouseChain';
            end
            obj.DataPreallocation % preallocate data structures           
            obj.SetFitOptions;    % set fit options
            
            obj.Run;
            obj.AnalyzeResults;
            obj.DisplayAnalyzedData
        end
        
        function SetFitOptions(obj)
            % Set the options for the fitting procedure
            obj.fitModel             = @(b,dists)(1/sum(dists.^(-b)))*dists.^(-b);
            obj.goalFunction         = @(b)sum((obj.fitModel(b,1:obj.params.numBeads-1)-obj.fitModel(1.5,1:obj.params.numBeads-1)).^2);
            obj.fOptions             = optimset;
            obj.fOptions.TolX        = 1e-16;
            obj.fOptions.TolFun      = 1e-16;
            obj.fOptions.MaxFunEvals = 2500;
            obj.fOptions.MaxIter     = 2500;
        end
        
        function DataPreallocation(obj)
            % preallocations
            obj.fitResults.mean           = cell(1,obj.params.numRounds);
            obj.fitResults.bead.fittedExp = zeros(obj.params.numBeads,obj.params.numRounds);
            obj.fitResults.bead.gof       = cell(obj.params.numBeads,obj.params.numRounds);
%             obj.beadDistance              = zeros(obj.params.numBeads,obj.params.numBeads,obj.params.numSimulations,obj.params.numRounds);
            obj.beadEncounterHistogram    = zeros(obj.params.numBeads,obj.params.numBeads-1,obj.params.numRounds);
            obj.encounterHistogram        = zeros(obj.params.numBeads,obj.params.numBeads,obj.params.numRounds);
            obj.beadEncounterHistogram    = zeros(obj.params.numBeads,obj.params.numBeads-1,obj.params.numRounds);
            obj.beadEncounterProbability  = zeros(obj.params.numBeads,obj.params.numBeads-1,obj.params.numRounds);
            obj.meanEncounterProbability  = zeros(1,obj.params.numBeads-1,obj.params.numRounds);
            obj.msd                       = zeros(obj.params.numBeads,obj.params.numSteps,obj.params.numRounds);
        end
               
        function ReadRecipeFile(obj)
            
            try
                % Read the content of the recipe file 
                t = fileread(fullfile(obj.params.recipeFolder, [obj.params.recipeFileName,'.rcp']));
            catch
                % If the recipe file was not found, switch to default
                % recipe
                warning('%s%s%s','Recipe file',[obj.params.recipeFileName,'.rcp'], ' was not found! Continuing with defaultRecipe') 
                t = fileread(fullfile(obj.params.recipeFolder, [obj.params.defaultRecipe,'.rcp']));
            end
            
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
        
        function SetRecipeParams(obj)
            % set the parameters defined in the recipe file 
            try % to support older versionof recipe files
             eval(obj.recipe.SetRecipeParams);
            catch
            end
        end
        
        function PreRoundActions(obj)
            obj.round = obj.round+1;
            eval(obj.recipe.PreSimulationBatchActions);
            
        end
        
        function PreRunActions(obj)
            % actions performed before each simulation run
            obj.simulation = obj.simulation+1;
            eval(obj.recipe.PreRunActions);            
            % initialize a new chain             
            if strcmpi(obj.polymerModel,'rouseChain')
                obj.handles.classes.chain = SimpleRouse(obj.params);           
            else 
                obj.handles.classes.chain = BetaPolymer(obj.params);
            end
              obj.handles.classes.chain.Initialize              
        end
        
        function PreStepActions(obj)
            % actions before each simulation step
            obj.step = obj.step+1;
            eval(obj.recipe.PreStepActions);
        end
        
        function Run(obj)
            % Run simulations                        
            for rIdx = 1:obj.params.numRounds
                t1Round = clock;
                obj.PreRoundActions;                
                for sIdx = 1:obj.params.numSimulations
                    t1Simulation= clock;
                    obj.PreRunActions                    
                    while obj.step<obj.params.numSteps && ~obj.stepExitFlag
                        obj.PreStepActions
                        obj.handles.classes.chain.Step;
                        obj.PostStepActions
                    end                                        
                    
                    % imediatly add to the encounter histogram 
                      obj.encounterHistogram(:,:,rIdx) = obj.encounterHistogram(:,:,rIdx)+double(obj.handles.classes.chain.beadDist(:,:,end)<obj.params.encounterDist);
                      obj.encounterHistogram(:,:,rIdx) = obj.encounterHistogram(:,:,rIdx)-diag(diag(obj.encounterHistogram(:,:,rIdx)));
                       
%                     obj.beadDistance(:,:,sIdx,rIdx) = obj.handles.classes.chain.beadDist(:,:,end);% take the last recorded bead distances
                    t2Simulation = clock;
                    obj.simulationData(obj.round).simulation(sIdx).time         = etime(t2Simulation,t1Simulation);
                    obj.simulationData(obj.round).simulation(sIdx).numSteps     = obj.handles.classes.chain.params.numSteps;
                    obj.simulationData(obj.round).simulation(sIdx).meanStepTime = obj.handles.classes.chain.simulationTime/obj.handles.classes.chain.params.numSteps;
                    obj.simulationData(obj.round).simulation(sIdx).params       = obj.params;
                    
                    obj.PostRunActions;% perform post run actions
                    fprintf('%s%d%s%d%s%s\n', 'round ', rIdx,' simulation ', sIdx, ' is done. Time ', num2str( obj.simulationData(obj.round).simulation(sIdx).time))
                end
                
                % record simulation data
                t2Round = clock;
                obj.simulationData(obj.round).round.time = etime(t2Round,t1Round);% record round time
                obj.simulationData(obj.round).params     = obj.params;% save round parameters
                
                % perform user defined actions at the end of each round
                obj.PostRoundActions
            end
        end
        
        function PostStepActions(obj)
            % actions performed after each step
            eval(obj.recipe.PostStepActions);
            obj.msd(:,obj.step,obj.round) = (obj.msd(:,obj.step,obj.round)*(obj.simulation-1) + obj.handles.classes.chain.msd)/obj.simulation;
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
            if obj.params.analyzeResults
%                 obj.CalculateEncounterHistogram;
                obj.CalculateBeadEncounterHistogram;
                obj.CalculateBeadEncounterProbability;
                obj.CalculateMeanEncounterProbability;
                obj.FitMeanEncounterProbability
                obj.FitBeadEncounterProbability;
                obj.FitMSD
            end            
        end
        
        function DisplayAnalyzedData(obj)
            if obj.params.analyzeResults
                obj.DisplayFittedExponents;
                obj.DisplayMeanEncounterProbability
                obj.DisplayBeadDataFit
                obj.DisplayEncounterHistograms
                obj.DisplayFittedAlphasOfMSD
%                 obj.DisplayMSD
            end
        end
        
        function CalculateEncounterHistogram(obj)%unused
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
%                 places = obj.meanEncounterProbability(1,:,rIdx)~=0;
                f =@(b)sum(((1/sum((dists).^(-b)))*(dists).^(-b)-obj.meanEncounterProbability(1,:,rIdx)').^2);
                % set zero weights to points with no data 
                [obj.fitResults.mean{rIdx}.beta] = fminbnd(f,0,2.2,obj.fOptions);
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
                    e = fminbnd(f,0,2.2,obj.fOptions);
                    
                    obj.fitResults.bead.fittedExp(bIdx,pIdx) = e;%beadFit.b;
                    %                 obj.fitResults.bead.fittedBias(bIdx,pIdx) = beadFit.a;
                    %                 obj.fitResults.bead.gof{bIdx,pIdx} = gof;
                end
            end
        end
        
        function FitMSD(obj)
            % estimate the exponent value for the msd of each bead                        
            if obj.params.calculateMSD
                f = fittype('A*x.^alpha');
%                 f = @(alpha) sum(
                opt = fitoptions(f);
                set(opt,'StartPoint',[0.5,0.5],...
                                 'Lower',[0 0]);
             for rIdx = 1:obj.params.numRounds
                for bIdx = 1:obj.params.numBeads
                    [fObj] = fit((1:numel(obj.msd(bIdx,:,rIdx)))',obj.msd(bIdx,:,rIdx)',f,opt);
                    obj.fitResults.bead.fittedAlphas(bIdx,rIdx)  = fObj.alpha;
                    obj.fitResults.bead.fittedAForMSD(bIdx,rIdx) = fObj.A;
                end
             end
            end
        end
        
        function DisplayFittedExponents(obj)
            % Display the values of beta obtained for each bead in each
            % round
            f = figure('Name','fittedExp','FileName','fittedExp');
            a = axes('Parent',f,'NextPlot','Add','FontSize',40);
            title(a,'Fitted \beta values','FontSize',40);
            for rIdx = 1:obj.params.numRounds           
             cColor = 0.9*rand(1,3);
             line('XData',1:obj.params.numBeads,...
                  'YData',obj.fitResults.bead.fittedExp(:,rIdx),...
                  'Marker','.',...
                  'MarkerSize',8,...                 
                  'LineWidth',7,...
                  'Color',cColor,...
                  'Parent',a,...
                  'DisplayName',['Experiment ', num2str(rIdx)]);
               % plot mean beta value
               
            line('XData',[1 obj.params.numBeads],...
                'YData',[obj.fitResults.mean{rIdx}.beta,obj.fitResults.mean{rIdx}.beta],...                
                'Color',cColor,...
                'Parent',a,...
                'LineWidth',4);
            end
           
            xlabel(a,'Bead','FontSize',40);
            ylabel(a,'\beta','FontSize',40);
            l = legend(flipud(get(a,'children')));
            set(l,'FontSize',30);
        end
        
        function DisplayBeadDataFit(obj,beadRange)
            % Display fit by bead for each round
            
            if ~exist('beadRange','var')
                beadRange = 1:obj.params.numBeads;
            end
                                    
            for rIdx = 1:obj.params.numRounds
                
                figName = sprintf('%s%s','BeadDataFitExperiment',num2str(rIdx));
                f       = figure('Name',figName,'FileName',figName);
                a(rIdx) = axes('NextPlot','Add','FontSize', 40,'Parent',f);% create axes                
                dists   = 1:obj.params.numBeads-1;
                for bIdx = 1:numel(beadRange)
                    lineColor = rand(1,3);
                    % plot the Encounter probability from simulation
                    line('Parent',a(rIdx),...
                        'XData',dists,...
                        'YData',obj.beadEncounterProbability(bIdx,:,rIdx),...
                        'Color',lineColor,...
                        'DisplayName',sprintf('%s%s','Bead',num2str(bIdx)),...
                        'LineWidth',5);
                    
                    % plot the fitted line
                    line('Parent',a(rIdx),...
                        'XData',dists,...
                        'YData',(1/sum(dists.^(-obj.fitResults.bead.fittedExp(bIdx,rIdx))))*dists.^(-obj.fitResults.bead.fittedExp(bIdx,rIdx)),...%obj.fitModel(obj.fitResults.bead.fittedBias(bIdx,rIdx),obj.fitResults.bead.fittedExp(bIdx,rIdx),dists),...
                        'Color',lineColor,...
                        'DisplayName',sprintf('%s%s','\beta= ',num2str(obj.fitResults.bead.fittedExp(bIdx,rIdx))),...
                        'LineWidth',5);
                end
                    % display expected model
                line('XData',dists,...
                     'YData',(1/sum(dists.^(-1.5*(mean(obj.params.beta)-1))))*dists.^(-1.5*(mean(obj.params.beta)-1)),...%obj.fitModel(-0.5/(-1+obj.params.numBeads^-0.5),1.5,dists),...
                     'Color','r',...
                     'Parent',a(rIdx),...
                     'DisplayName',['Theoretical encounter \beta=' num2str(1.5*(mean(obj.params.beta)-1))],...
                     'LineWidth',5)
                
                title(a(rIdx),sprintf('%s%s','BeadDataFitExperiment',num2str(rIdx)),'FontSize',40)
                xlabel(a(rIdx),'Bead distance','FontSize',40);
                ylabel(a(rIdx),'Encounter prob','Fontsize',40);
            end
        end
        
        function DisplayMeanEncounterProbability(obj)
            f = figure('Name','meanEncounterProbability',...
                       'FileName','meanEncounterProbability') ;
            a = axes('Parent',f,'NextPlot','Add','FontSize',40);
            
            title(a,'Mean Encounter Prob.','FontSize',40);
            xlabel(a,'Distance [bead]','FontSize',40);
            ylabel(a,'Encounter Prob.','FontSize',40);
            
            dists = (1:obj.params.numBeads-1)';
            for rIdx = 1:obj.params.numRounds
                c = rand(1,3);
                line('XData',dists,...
                    'YData',obj.meanEncounterProbability(1,:,rIdx),...
                    'Color',c,...
                    'Parent',a,...
                    'LineWidth',5,...
                    'DisplayName',['Round ', num2str(rIdx)])
                
                line('XData',dists,...
                    'YData',(1/sum(dists.^(-obj.fitResults.mean{rIdx}.beta)))*dists.^(-obj.fitResults.mean{rIdx}.beta),...%obj.fitModel(obj.fitResults.mean{rIdx}.a,obj.fitResults.mean{rIdx}.b,dists),...
                    'Color',c,...
                    'DisplayName',[ '\beta=', num2str(obj.fitResults.mean{rIdx}.beta)],...
                    'Parent',a,...
                    'LineWidth',5)
                %                 m(rIdx) = obj.fitResults.mean{rIdx}.a;
            end
            % add the theoretical model (temporary)
            
            line('XData',dists,...
                'YData',(1/sum(dists.^(-1.5*(mean(obj.params.beta)-1))))*dists.^(-1.5*(mean(obj.params.beta)-1)),...
                'Color','r',...
                'Parent',a,...
                'LineWidth',5,...
                'DisplayName',['Theoretical encounter \beta=', num2str(1.5*(mean(obj.params.beta)-1))])
           l= legend(flipud(get(a,'Children')));
           set(l,'FontSize',18);

        end
        
        function DisplayEncounterHistograms(obj)
            
            for rIdx = 1:obj.params.numRounds
                figure('Name',['Experiment',num2str(rIdx)],'FileName',['EncounterHistogramExperiment',num2str(rIdx)]);
                windowSize = 10;
                imagesc(medfilt2(obj.encounterHistogram(:,:,rIdx),[windowSize,windowSize])), colormap hot
                title(sprintf('%s%d','Experiment ', rIdx),'FontSize',40); 
                xlabel('Bead number','FontSize',40);
                ylabel('BeadNumber','FontSize', 40);
                set(gca,'FontSize',40);
                daspect([1 1 1])
            end
            
        end
        
        function DisplayFittedAlphasOfMSD(obj)
            
            for rIdx = 1:obj.params.numRounds
                f = figure('Units','norm',...
                    'Name',['\alpha values experiment' num2str(rIdx)],...
                    'FileName',['AlphaValuesOfMSDExperiment' num2str(rIdx)]);
                a= axes('Parent',f,'FontSize',40,'LineWidth',4,'NextPlot','Add');
                line('XData',1:obj.params.numBeads,...
                    'YData',obj.fitResults.bead.fittedAlphas(:,rIdx),...
                    'Marker','o',...
                    'MarkerSize',7,...
                    'MarkerFaceColor','b',...
                    'MarkerEdgeColor',[0.5 0.1 0.2],...
                    'LineWidth',4,...
                    'Parent',a);
                title(a,'Fitted alpha','FontSize',40);
                xlabel(a,'bead','FontSize',40),
                ylabel(a,'\alpha','FontSize',40,'Rotation',0);   
                
                % Add the mean line 
                m = mean(obj.fitResults.bead.fittedAlphas(:,rIdx));
                line('XData',[1 obj.params.numBeads],...
                     'YData',[m m],...
                     'LineWidth',4,...
                     'LineStyle','-.',...
                     'Color','g',...
                     'Parent',a,...
                     'DisplayName','mean \alpha')
                
            end
        end
        
        function DisplayMSD(obj)
            if obj.params.calculateMSD
             for rIdx = 1:obj.params.numRounds
                 f = figure('Units','norm','Name','MSD','FileName',['MsdExperiment', num2str(rIdx)]);
                 a = axes('Parent',f,'fontSize',40,...
                     'NextPlot','Add');
                 t = 1:numel(obj.msd(1,:,1));
                 for bIdx = 1:obj.params.numBeads
                  c = rand(1,3);
                  line('XData',t,...
                      'YData',obj.msd(bIdx,:,rIdx),...
                      'LineWidth',3,...
                      'Color',c,...
                      'Parent',a,...
                      'DisplayName',['bead ' num2str(bIdx)])
                  line('XData',t,...
                       'YData',obj.fitResults.bead.fittedAForMSD(bIdx,rIdx)*t.^obj.fitResults.bead.fittedAlphas(bIdx,rIdx),...
                       'Color','r',...
                       'DisplayName',['\alpha = ' num2str(obj.fitResults.bead.fittedAlphas(bIdx,rIdx))])
                 end            
             end
            end
        end
        
    end
end