classdef SimulationDataRecorder<handle
    % Holds data coming from Rouse chain simulations. This class replaces
    % the DistributionHandler in collecting the data
    properties
        simulationInfo       % describes the simulation
        simulationData  = struct('chainObj',[],...
            'numChains',0,...
            'step',0,...
            'time',0,...
            'positions',[],...
            'beadDist',cell(1)); % stores data regarding each simulation step
        simulationBatchRound = 0  % the number of simulation batch
        simulationRound      = 0; % the simulation index
        params
    end
    
    events
    end
    
    methods
        function obj = SimulationDataRecorder(dataRecorderParams)
            % class constructor
            obj.LoadDefaultParams
            obj.LoadInputParams(dataRecorderParams);
        end
        
        function LoadDefaultParams(obj)
            obj.params.saveType = 'external';% save data to external mat files. no data is saved on the class
            obj.params.resultsFolder = pwd; % default results folder
        end
        
        function LoadInputParams(obj,params)
            f = fieldnames(params);
            for fIdx =1:numel(f)
                obj.params.(f{fIdx}) = params.(f{fIdx});
            end
        end
        
        function Add(obj)
            % Add data to the stored distributions
            % each call to Add raises the step count
            obj.simulationData(obj.simulationRound).step = obj.simulationData(obj.simulationRound).step+1;
            obj.AddBeadDist
            obj.AddBeadsPosition
            obj.AddBeadEncounterHist
        end
        
        function AddBeadDist(obj)
            % Add the bead distance matrix to the current simulation data
            for cIdx = 1:obj.simulationData(obj.simulationRound).numChains
                obj.simulationData(obj.simulationRound).beadDist{cIdx}(:,:,obj.simulationData(obj.simulationRound).step) = ...                   
                    obj.simulationData(obj.simulationRound).chainObj(cIdx).beadsDist;
                obj.simulationData(obj.simulationRound).beadDistSquare{cIdx}(:,:,obj.simulationData(obj.simulationRound).step) =...
                    obj.simulationData(obj.simulationRound).chainObj(cIdx).beadsDist.^2;
            end
        end
        
        function AddBeadEncounterHist(obj)
            % find beads which are closer than the encounter distance
            % defined        
            for cIdx = 1:obj.simulationData(obj.simulationRound).numChains
                f = double(obj.simulationData(obj.simulationRound).beadDist{cIdx}(:,:,obj.simulationData(obj.simulationRound).step)<obj.params.encounterDist);   
                f= f-eye(size(f));
                % add simulation first encounter time 
%                 eTime = double(f).*obj.simulationData(obj.simulationRound).step;
%                 eTime(eTime==0) = NaN;                
%                 obj.simulationData(obj.simulationRound).encounterTime{cIdx}(:,:,obj.simulationData(obj.simulationRound).step) = eTime;
                obj.simulationData(obj.simulationRound).encounterHist{cIdx} = obj.simulationData(obj.simulationRound).encounterHist{cIdx}+f;
            end
        end
        
        function AddBeadsPosition(obj)
            % Record the position of the rouse chains
            for cIdx = 1:obj.simulationData(obj.simulationRound).numChains
                obj.simulationData(obj.simulationRound).positions(cIdx) = ...
                    obj.simulationData(obj.simulationRound).chainObj(cIdx).positions.beads.cur;
            end
        end
        
        function NewSimulation(obj,chainObj,frameworkParams)
            obj.CreateNewSimulationStruct(chainObj,frameworkParams);
        end
        
        function CreateNewSimulationBatch(obj)
            % Activated upon new simulation batch
            obj.simulationBatchRound = obj.simulationBatchRound+1;% advance the count
            obj.simulationRound      = 0;% reset the simulation round counter;
        end
        
        function CreateNewSimulationStruct(obj,chainObj,frameworkParams)
            % Initialize new simulation structure witht he chain object
            % participating in the simulation round
            
            obj.simulationRound = obj.simulationRound+1; % advance the index by 1
            
            obj.simulationData(obj.simulationRound).chainObj  =  chainObj;             
            obj.simulationData(obj.simulationRound).numChains =  numel(chainObj);
            obj.simulationData(obj.simulationRound).step      =  1;
            obj.simulationData(obj.simulationRound).time      =  0;
            obj.simulationData(obj.simulationRound).parameters = frameworkParams;
            for cIdx = 1:obj.simulationData(obj.simulationRound).numChains
                % add initialized beadDistance matrix
                obj.simulationData(obj.simulationRound).beadDist{cIdx} = chainObj(cIdx).beadsDist;
                obj.simulationData(obj.simulationRound).encounterHist{cIdx} = zeros(chainObj(cIdx).params.numBeads);
                obj.simulationData(obj.simulationRound).encounterTime{cIdx} = zeros(chainObj(cIdx).params.numBeads);
                % record positions
                obj.simulationData(obj.simulationRound).positions.x = chainObj(cIdx).positions.beads.cur.x;
                obj.simulationData(obj.simulationRound).positions.y = chainObj(cIdx).positions.beads.cur.y;
                obj.simulationData(obj.simulationRound).positions.z = chainObj(cIdx).positions.beads.cur.z;
            end
        end
        
        function SetSimulationStartTime(obj)
            % Record the time the simulation started
            startTime = clock;
            obj.simulationData(obj.simulationRound).startTime.year   = startTime(1);
            obj.simulationData(obj.simulationRound).startTime.month  = startTime(2);
            obj.simulationData(obj.simulationRound).startTime.day    = startTime(3);
            obj.simulationData(obj.simulationRound).startTime.hour   = startTime(4);
            obj.simulationData(obj.simulationRound).startTime.minute = startTime(5);
            obj.simulationData(obj.simulationRound).startTime.second = startTime(6);
        end
        
        function SetSimulationEndTime(obj)
            % Record the time the simulation ended
            endTime = clock;
            obj.simulationData(obj.simulationRound).endTime.year   = endTime(1);
            obj.simulationData(obj.simulationRound).endTime.month  = endTime(2);
            obj.simulationData(obj.simulationRound).endTime.day    = endTime(3);
            obj.simulationData(obj.simulationRound).endTime.hour   = endTime(4);
            obj.simulationData(obj.simulationRound).endTime.minute = endTime(5);
            obj.simulationData(obj.simulationRound).endTime.second = endTime(6);
            
            f = fieldnames(obj.simulationData(obj.simulationRound).startTime);
            startTime = zeros(1,numel(f));
            for fIdx = 1:numel(f)
                startTime(fIdx) = obj.simulationData(obj.simulationRound).startTime.(f{fIdx});
            end
            obj.simulationData(obj.simulationRound).simulationTime = etime(endTime,startTime);
        end
        
        function SaveResults(obj)
            % Save results. implement 4 type of saving
            % [internal/external/none/all].
            if strcmpi(obj.params.saveType,'external')
                % Save to external mat files, no data is saved in the
                % class. use this option when memery problems occur
                obj.ExportData
%                 obj.ClearCurrentSimulationData   % clear data from fields
            elseif strcmpi(obj.params.saveType,'none')
                % Clear data from fields dont export data to mat files
                obj.ClearCurrentSimulationData;
            elseif strcmpi(obj.params.saveType,'all')
                % Export data to mat file while keeping the data internally.
                % Beware of memory problems for long simulations
                obj.ExportData
            elseif strcmpi(obj.params.saveType,'internal')
                % Don't export data to mat files. keep the data in the class properties
            end
        end
        
        function ClearCurrentSimulationData(obj)
            % Clear the data from the fields
            for cIdx = 1:obj.simulationData(obj.simulationRound).numChains
                obj.simulationData(obj.simulationRound).positions(end).x = [];
                obj.simulationData(obj.simulationRound).positions(end).y = [];
                obj.simulationData(obj.simulationRound).positions(end).z = [];
                obj.simulationData(obj.simulationRound).chainObj(end)    = [];
                obj.simulationData(obj.simulationRound).numChains        = [];
                obj.simulationData(obj.simulationRound).step             = [];
                obj.simulationData(obj.simulationRound).time             = [];
                obj.simulationData(obj.simulationRound).beadDist{end}    = [];
            end
        end
        
        function ClearAllSimulationData(obj)
            for cIdx = 1:obj.simulationData(obj.simulationRound).numChains
                for rIdx = 1:obj.simulationRound
                    obj.simulationData(rIdx).positions(cIdx).x = [];
                    obj.simulationData(rIdx).positions(cIdx).y = [];
                    obj.simulationData(rIdx).positions(cIdx).z = [];
                    obj.simulationData(rIdx).chainObj(cIdx)    = [];
                    obj.simulationData(rIdx).numChains         = [];
                    obj.simulationData(rIdx).step              = [];
                    obj.simulationData(rIdx).time              = [];
                    obj.simulationData(rIdx).beadDist{cIdx}    = [];
                    obj.simulationData(rIdx).simulationTime    = 0;
                end
            end
%                         obj.simulationData = struct('chainObj',[],...
%                                     'numChains',0,...
%                                     'step',0,...
%                                     'time',0,...
%                                     'positions',[],...
%                                     'beadDist',cell(1)); % stores data regarding each simulation step
        end
        
        function ExportData(obj)
            % Export data to external mat file
            c = clock;
            currentSimulationResultsPath = sprintf('%s%s%s_%s_%s',obj.params.resultsFolder,'\',num2str(c(3)), num2str(c(2)), num2str(c(1)));
            [~,~,~]  = mkdir(currentSimulationResultsPath);
            for cIdx = 1:obj.simulationData(obj.simulationRound).numChains
                fileName = sprintf('%s%s%s%s%s%s','SimulationBatch_',num2str(obj.simulationBatchRound),...
                                                  '_SimulationRound',num2str(obj.simulationRound),...
                                                  '_Chain_',num2str(cIdx));
                results.beadDistMat       = obj.simulationData(obj.simulationRound).beadDist{cIdx}(:,:,3000:end); % save the mean distance matrix
                results.beadDistanceRMS   = sqrt(mean(obj.simulationData(obj.simulationRound).beadDistSquare{cIdx}(:,:,3000:end),3));
                results.beadEncounterHist = obj.simulationData(obj.simulationRound).encounterHist{cIdx};
                results.beadEncounterTime = obj.simulationData(obj.simulationRound).encounterTime{cIdx};
                results.params            = obj.simulationData(obj.simulationRound).parameters;
                save(fullfile(currentSimulationResultsPath,fileName),'results','-v7.3');
            end
            
            % save the parameter file and the recipe file 
            copyfile(fullfile(pwd,'Framework','SimulationFrameworkParams.xml'),fullfile(currentSimulationResultsPath,'SimulationFrameworkParams.xml'));
            copyfile(fullfile(pwd,'Framework','Recipes',[obj.params.recipeFileName '.rcp']),currentSimulationResultsPath);
        end
    end
end
