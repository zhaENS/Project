classdef RouseSimulatorFramework<handle      
    %TODO: add presets for simulations,(load/save)
    %TODO: move the initialization of the chains to allow initialization on
    %      the fly
    %TODO: add a listener that draws emails from the server, to allow
    %      controlling the simulations from a distance
    %TODO: save noise terms seeds to allow re-simulations
    %TODO: save parameters and a short explanation with the results   
    %TODO: remove the option to plot online. plotting is always offline.
    %      remove parentAxes property from domainHandlerParams
    %TODO: chain parameters should be initialized before the chains according to the number specified by the simulator params 
    %TODO: build a simple GUI
    %TODO: insert the correct indices for fixed particles in Step 
    %TODO: fix reflection for general meshes shapes
    %TODO: consider moving Domain to the ObjectManager
    % TODO: find a uniform framework for domains and chains
    % TODO: integrate the domain with the object manager.
    % TODO: allow connection of sticky beads
    % TODO: expand ChainParams.m to include parameters for attaching and detaching either from domain or from other chains
    % TODO: complete diffusion constrained on domain
    % TODO: complete diffusion constrained on chains
    % TODO: add documentation
    % TODO: check graph representation for object mapper
    % TODO: allow some bonds to be stiff, find a way to simulate a stiff connector chain efficiently
    % TODO: add reflect outside the domain
    % TODO: add constrained diffusing objects to exert force on the objects
    properties
        handles     
        objectManager      % data and object manipulation
        dataRecorder       % simulation data module
        simulationGraphics % graphical module
        runSimulation      % flag 
        params             % object parameters
        simulationData  = struct('step',[],'stepTime',[]); % simulation data (should move to dataRecorder)
        batchRound      = 0;
        simulationRound = 0; 
    end
    
    properties (Access=private)
        chainsCenterOfMass
        recipe = struct('PreSimulationBatchActions','',...
                        'PreRunActions','',...
                        'PostRunActions','',...
                        'PostSimulationBatchActions','');
                    
       SetRecipeParamsFlag            = false;
       PreSimulationBatchActionsFlag  = false;
       PreRunActionsFlag              = false;
       PostRunActionsFlag             = false;
       PostSimulationBatchActionsFlag = false;
       
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
            obj.params = frameworkParams;           

        end
        
        function InitializeClasses(obj)
            
            obj.InitializeDomain;                        
            
            obj.InitializeObjectManager;
                        
            obj.InitializeDataRecorder;
            
            obj.InitializeGraphics;
        end                              
        
        function InitializeDataRecorder(obj)
            obj.dataRecorder = SimulationDataRecorder(obj.params.dataRecorder);
        end
        
        function InitializeGraphics(obj)
          % Initialize framework's graphics 
          obj.simulationGraphics = SimulationFrameworkGraphics(obj);   
          obj.simulationGraphics.CreateControls;
        end
        
        function InitializeDomain(obj)
             % Parameters for the domain are set in the organizeParams             
             obj.handles.classes.domain  = DomainHandler(obj.params.domain);
        end   
        
        function InitializeObjectManager(obj)
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
                        % advance one step of the polymer chain    
                        obj.Step;                       
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
            obj.dataRecorder.CreateNewSimulationBatch
            
            % evaluate the recipe file at the PreSimulationBatchActions
            eval(obj.recipe.PreSimulationBatchActions)
        end
        
        function PreRunActions(obj)
            % Actions called before running the simulation        
            obj.simulationRound = obj.simulationRound+1;
            obj.runSimulation   = true;            
            
            eval(obj.recipe.PreRunActions);
            
            obj.dataRecorder.NewSimulation(obj.objectManager.handles.chain,obj.params);%(obj.handles.classes.rouse,obj.params);                        
            obj.dataRecorder.SetSimulationStartTime;
            obj.simulationData(obj.batchRound,obj.simulationRound).step     = 1;
            obj.simulationData(obj.batchRound,obj.simulationRound).stepTime = 0;
        end
        
        function PreStepActions(obj)
            eval(obj.recipe.PreStepActions);
        end
        
        function Step(obj,varargin)%TODO: consider joining Domain to objectManager
            % Next simulation step 
            objList        = 1:obj.objectManager.numObjects;
                        
            % Advance one step and apply object forces                        
            obj.objectManager.Step(objList)% update current and previous object position
            
            % Update object list
            objList        = 1:obj.objectManager.numObjects;
            
            % Apply domain (global) forces on all objects in the domain 
            dp                      = obj.handles.classes.domain.params.forceParams;% Get domain parameters
            prevParticlePosition    = obj.objectManager.prevPos;      % prev position 
            curParticlePosition     = obj.objectManager.curPos;       % new pos after internal forces
            particleDist            = obj.objectManager.particleDist; % distance before applying internal forces
            fixedParticleNum        = [obj.objectManager.fixedParticles{:}];% ToDo: insert the correct indices
            
            % Apply external forces from the domain and reflect
            curParticlePosition = obj.handles.classes.domain.ApplylForces(prevParticlePosition,...
                                                                          curParticlePosition,...
                                                                          particleDist,...
                                                                          dp.lennardJonesForce,dp.diffusionForce,...
                                                                          dp.diffusionConst,...    
                                                                          dp.LJPotentialWidth,dp.LJPotentialDepth,...
                                                                          fixedParticleNum,obj.params.simulator.dt);                                         
                                                       
            % Deal the positions after reflection between the objects and their members in
            % the domain 
            obj.objectManager.DealCurrentPosition(objList,curParticlePosition);
            obj.objectManager.DealPreviousPosition(objList,curParticlePosition);
                                    
            % Show simulation
            obj.simulationGraphics.ShowSimulation
            
            % check for interaction between objects             
%              obj.objectManager.ObjectInteraction;
             
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
            
            obj.dataRecorder.SetSimulationEndTime;
            obj.dataRecorder.SaveResults;% save results
            
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
            obj.dataRecorder.ClearAllSimulationData;            
%             obj.simulationRound = 0;               
        end                                        
        
        function Record(obj)%TODO: verify data passes correctly 
            % Add samples to the recorded distributions in
            % distributionHandler
            if obj.params.simulator.recordData
               obj.dataRecorder.Add;
            end
        end                                                
        
    end
       
end
    