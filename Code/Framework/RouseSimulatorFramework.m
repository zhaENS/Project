classdef RouseSimulatorFramework<handle      
    %TODO: add presets for simulations,(load/save)
    %TODO: move the initialization of the chains to allow initialization on
    %      the fly
    %TODO: add a listener that draws emails from the server, to allow
    %controlling the simulations from a distance
    %TODO: save noise terms seeds to allow resimulations
    %TODO: save parameters and a short explanation with the results   
    %TODO: remove the option to plot online. plotting is always offline.
    %      remove parentAxes property from domainHandlerParams
    %TODO: chain parameters should be initialized before the chains according to the number specified by the simulator params     
    properties
        handles     
        objectManager
        beadDistance
        runSimulation
        params
        simulationData  = struct('step',[],'stepTime',[]);      
        batchRound      = 0;
        simulationRound = 0; 
    end
    
    properties (Access=private)
%         chainColors 
        chainsCenterOfMass
        recipe = struct('PreSimulationBatchActions','',...
                        'PreRunActions','',...
                        'PostRunActions','',...
                        'PostSimulationBatchActions','');
    end    
    
    methods 
        function obj = RouseSimulatorFramework(frameworkParams)
                                
           % Set initial parameters
           obj.OrganizeParams(frameworkParams)

           % Initizlie classes
           obj.InitializeClasses
           
           % set recipe files 
           obj.ReadRecipeFile           
        end
        
        function OrganizeParams(obj,frameworkParams)% move to params class
            obj.params                             = frameworkParams;           
%             obj.params.chain.dt                    = obj.params.simulator.dt;
%             obj.params.chain.dimension             = obj.params.simulator.dimension;
%             obj.params.domain.showDomain           = obj.params.simulator.showSimulation;  
%             obj.params.domain.dimension            = obj.params.simulator.dimension;
%             obj.params.domain.dt                    = obj.params.simulator.dt;
%             obj.params.dataRecorder.recipeFileName = obj.params.simulator.recipeFileName;
%             obj.params.dataRecorder.encounterDist  = obj.params.simulator.encounterDist;
        end
        
        function InitializeClasses(obj)
            
            obj.InitializeDomain;                        
            
            obj.InitializeObjects;
                        
            obj.InitializeDataRecorder;
            
            obj.InitializeGraphics;
        end
        
        function SetInputParams(obj,varargin)
            if mod(numel(varargin),2)~=0
                error('value pair input must come in pairs')
            end
            for pIdx = 1:2:numel(varargin);
                obj.params.(varargin{pIdx}) = varargin{pIdx+1};
            end
        end                        
        
        function InitializeGraphics(obj)
          % Initialize framework's graphics 
          obj.handles.classes.graphics = SimulationFrameworkGraphics(obj);   
          obj.handles.classes.graphics.CreateControls;
        end
        
        function InitializeObjects(obj)
            % Initialize objectManager
            obj.objectManager = ObjectManager(obj.params.chain);
            obj.objectManager.InitializeObjects(obj.handles.classes.domain);
                                                                               
        end        
        
        function ReadRecipeFile(obj)
            
         t = fileread([obj.params.simulator.recipeFileName,'.rcp']);
         % Search for the function marker 
         [funcStartPos1,funcStartPos2] = regexp(t,'<func>');
         [funcEndPos1,funcEndPos2]     = regexp(t,'</func>');
         for fIdx = 1:numel(funcStartPos1)-1
             % Sort functions into categories
             functionName = strtrim(t(funcStartPos2(fIdx)+1:funcEndPos1(fIdx)-1));
             obj.recipe.(functionName) = t(funcEndPos2(fIdx)+1:funcStartPos1(fIdx+1)-1);
         end
             functionName = strtrim(t(funcStartPos2(end)+1:funcEndPos1(end)-1));
             obj.recipe.(functionName) = t(funcEndPos2(end)+1:end);     
        end
               
        function InitializeDataRecorder(obj)
            obj.handles.classes.dataRecorder = SimulationDataRecorder(obj.params.dataRecorder);
        end
        
        function Run(obj)
            % Run the simulation according to the specified simulationType
            % parameter
           
            for bIdx = 1:obj.params.simulator.numSimulationBatches
                % perform actions before each simulation batch
                obj.PreSimulationBatchActions
                for sIdx = 1:obj.params.simulator.numSimulations
                    % perform actions before each simulation
                    obj.PreRunActions                    
                    while obj.runSimulation
                        % perform actions before the current step
                        obj.PreStepActions                                               
                        obj.Step;% advance one step of the polymer chain                           
                        % perform action post the current step 
                        obj.PostStepActions                                              
                    end
                    % perform actions post simulation
                    obj.PostRunActions
                end
                % perform actions post simulation batch
                obj.PostSimulationBatchActions
            end
           
        end         
        
        function PreSimulationBatchActions(obj)
            % Actios performed before a simulation batch      
            obj.batchRound      = obj.batchRound+1;
            obj.simulationRound = 0;
            obj.handles.classes.dataRecorder.CreateNewSimulationBatch
            
            % evaluate the recipe file at the PreSimulationBatchActions
            eval(obj.recipe.PreSimulationBatchActions)
        end
        
        function PreRunActions(obj)
            % Actions called before running the simulation        
            obj.simulationRound = obj.simulationRound+1;
            obj.runSimulation   = true;            
            
            eval(obj.recipe.PreRunActions);
            
            obj.handles.classes.dataRecorder.NewSimulation(obj.objectManager.handles.chain,obj.params);%(obj.handles.classes.rouse,obj.params);                        
            obj.handles.classes.dataRecorder.SetSimulationStartTime;
            obj.simulationData(obj.batchRound,obj.simulationRound).step     = 1;
            obj.simulationData(obj.batchRound,obj.simulationRound).stepTime = 0;
        end
        
        function PreStepActions(obj)
            eval(obj.recipe.PreStepActions);
        end
        
        function Step(obj,varargin)
            % Next simulation step 
            numObjects        = obj.objectManager.numObjects;
                        
            % Apply object forces
            for oIdx = 1:numObjects                                
                obj.objectManager.Step(oIdx)% update current object position
            end
            
            % Apply domain (global) forces            
            dp                = obj.handles.classes.domain.params.forceParams;% Get domain parameters
            [prevParticlePosition,curParticlePosition,connectivityMap,cp]= obj.objectManager.GetObjectsAsOne(1:numObjects);
            
            curParticlePosition = obj.handles.classes.domain.ApplyForces(curParticlePosition,connectivityMap,...                
                                             cp.springConst,dp.diffusionConst,cp.bendingConst,...
                                             dp.LJPotentialWidth,dp.LJPotentialDepth,...
                                             cp.minBeadDistance,cp.fixedBeadNum,obj.params.simulator.dt);
             % Reflect particles 
             [~,newPos] = obj.handles.classes.domain.Reflect(prevParticlePosition,...
                                                           curParticlePosition);
                                                       
            % Deal the positions after reflection between the objects in
            % the domain 
            obj.objectManager.DealCurrentPosition(1:numObjects,newPos);
            obj.objectManager.DealPreviousPosition(1:numObjects,newPos);
            prob = 0.995;
            for oIdx = 1:numObjects
                r = rand(1);
                if r>prob
                    obj.objectManager.ConnectParticles(oIdx,1,64);
                elseif r<(1-prob)
                    obj.objectManager.DisconnectParticles(oIdx,1,64)
                end
            end
%             for oIdx = 1:numObjects    
%                 % Get object data
%                 [prevPos, curPos]  = obj.objectManager.GetPosition(oIdx);
%                  cMap              = obj.objectManager.GetConnectionMap(oIdx);
%                  objParams         = obj.objectManager.GetObjectParameters(oIdx);
%                  
%                                  
%                 prevParticlePosition = prevPos{1};
%                 curParticlePosition  = curPos{1};
%                 
%                 connectivityMap     = cMap{1}.map;
%                 cp                  = objParams{1};                                          
%                 curParticlePosition = obj.handles.classes.domain.ApplyForces(curParticlePosition,connectivityMap,...                
%                                              cp.springConst,dp.diffusionConst,cp.bendingConst,...
%                                              dp.LJPotentialWidth,dp.LJPotentialDepth,...
%                                              cp.minBeadDistance,cp.fixedBeadNum,obj.params.simulator.dt);
%                                       
%                 
%                 % Reflect particles 
%                 [~,newPos]= obj.handles.classes.domain.Reflect(prevParticlePosition,...
%                                                            curParticlePosition);
%                                                        
%                 obj.objectManager.SetCurrentParticlePosition(oIdx,newPos); 
%                 obj.objectManager.SetPreviousParticlePosition(oIdx,newPos); 
%             end
             
            % Show simulation
            obj.handles.classes.graphics.ShowSimulation
            
            % Update simulation data
            obj.simulationData(obj.batchRound,obj.simulationRound).step = ...
                obj.simulationData(obj.batchRound,obj.simulationRound).step+1;
        end  
                
        function PostStepActions(obj)
            eval(obj.recipe.PostStepActions);
              obj.Record % record data related to the Rouse polymer                        
              % raise the Stop flag if the number of steps is more than the allowed
              obj.runSimulation = obj.runSimulation && ...
             (obj.simulationData(obj.batchRound,obj.simulationRound).step<obj.params.simulator.numSteps);
        end
        
        function PostRunActions(obj)
            % Calculate distributions related to the Rouse polymer
            % chain
            
            eval(obj.recipe.PostRunActions);
            
            obj.handles.classes.dataRecorder.SetSimulationEndTime;
            obj.handles.classes.dataRecorder.SaveResults;% save results
            
            if obj.params.simulator.notifyByEmail
                if mod(obj.simulationRound,obj.params.simulator.notifyCycleLength)==0
                    msg = sprintf('%s%s%s\n%s%s%s\n%s%s\n%s','simulation batch',num2str(obj.batchRound),', ',...
                        'simulation ', num2str(obj.simulationRound),' is done. ',...
                        'Simulation Time: ',num2str(obj.handles.classes.dataRecorder.simulationData(obj.simulationRound).simulationTime),...
                        'Simulating the distance between monomenrs at relexation time, each time bead 1 is connected to bead j, where j=3..64');
                    try
                        SendMail('ofirshukron','Simulation Status',msg);
                    catch
                    end
                end
            end
        end
        
        function PostSimulationBatchActions(obj)
            % actions performed post simulation batch
            eval(obj.recipe.PostSimulationBatchActions);            
            obj.handles.classes.dataRecorder.ClearAllSimulationData;            
%             obj.simulationRound = 0;               
        end        
        
        function InitializeDomain(obj)
            
             % parameters for the domain are set in the organizeParams
             
             obj.handles.classes.domain  = DomainHandler(obj.params.domain);
             
%             if obj.params.domain.showDomain
%                 % move these commands to future Plotter class               
%                obj.params.domain.parentAxes     = obj.handles.graphical.mainAxes;
%                obj.handles.classes.domain.ConstructDomain;
%                set(obj.handles.graphical.mainAxes,'NextPlot','Add')                 
%                set(obj.handles.graphical.mainAxes,'NextPlot','ReplaceChildren')
%             end
        end                           
        
        function Record(obj)%TODO: verify data passage 
            % Add samples to the recorded distributions in
            % distributionHandler
            if obj.params.simulator.recordData
%                obj.handles.classes.distributionHandler.Add% add to the frequency array of the distributionHandler
               obj.handles.classes.dataRecorder.Add;
            end
        end                                                
        
    end
       
end
    