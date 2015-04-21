classdef SimulationFrameworkParams<handle
    % this class holds all parameters used by the Framework class
    properties
        simulator
        chain        = ChainParams;
        domain       = DomainHandlerParams      
        dataRecorder = SimulationDataRecorderParams;
        plotHandler
    end
    
    methods
        
        function obj = SimulationFrameworkParams
            % class constructor 
            obj.SetSimulatorParams;
            obj.SetChainParams;
            obj.SetDomainParams;
            obj.SetDataRecorderParams;
            obj.SetPlotHandlerParams;
        end
        
        function SetSimulatorParams(obj)
            % Simulator framework params
            obj.simulator.dimension            = 3;
            obj.simulator.runSimulation        = false; % a flag indicating whether to allow the simulation to run at initiation 
            obj.simulator.numSimulationBatches = 1;     % number of simulation batches
            obj.simulator.numSimulations       = 1;     % number of simulations in each batch
            obj.simulator.numSteps             = 1000;   % for inf place inf
            obj.simulator.dt                   = 1e-2;  % time step 
            obj.simulator.numChains            = 5;   
            obj.simulator.encounterDist        = 0.1;   % The distance for which two monomer are considered to have met 
            obj.simulator.showSimulation       = true; 
            obj.simulator.recordData           = false;
            obj.simulator.notifyByEmail        = false;
            obj.simulator.notifyCycleLength    = 32;    % number of simulation cycles after which an email is sent 
            obj.simulator.recipeFileName       = 'debugRecipe';
            obj.simulator.recipesFolder        = ''; 
            
            % Control domain forces
            obj.simulator.diffusionConst       = 0.2;
            obj.simulator.LJPotentialWidth     = 0.01;
            obj.simulator.LJPotentialDepth     = 0.01;
        end
        
        function SetChainParams(obj)%TODO: change force params accordingly
            % Assign parameters for each chain 
            for cIdx = 1:obj.simulator.numChains
                obj.chain(cIdx) = ChainParams;
                % inherit the framework parameters 
                obj.chain(cIdx).diffusionConst = obj.simulator.diffusionConst;
                obj.chain(cIdx).dimension = obj.simulator.dimension;            
                obj.chain(cIdx).dt        = obj.simulator.dt;
                obj.chain(cIdx).SetForceParams;
            end
        end
        
        function SetDomainParams(obj)% TODO: change force params accordingly 
            obj.domain = DomainHandlerParams;
            
            % inherit the framework parameters 
            obj.domain.dimension        = obj.simulator.dimension;  
            obj.domain.diffusionConst   = obj.simulator.diffusionConst;
            obj.domain.showDomain       = obj.simulator.showSimulation;  
            obj.domain.dt               = obj.simulator.dt;
            obj.domain.LJPotentialWidth = obj.simulator.LJPotentialWidth;
            obj.domain.LJPotentialDepth = obj.simulator.LJPotentialDepth;
            obj.domain.SetForceParams;
           
        end
        
        function SetDataRecorderParams(obj)
            obj.dataRecorder = SimulationDataRecorderParams;
            % inherit the framework params
            obj.dataRecorder.recipeFileName = obj.simulator.recipeFileName;
            obj.dataRecorder.encounterDist  = obj.simulator.encounterDist;
        end
        
        function SetPlotHandlerParams(obj)
        end
        
    end
end