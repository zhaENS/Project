classdef SimulationFrameworkParams<handle
    % this class holds all parameters used by the Framework class
    properties
        simulator
        chain
        domain       
        dataRecorder
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
            obj.simulator.numSteps             = 500;   % for inf place inf
            obj.simulator.numChains            = 1;   
            obj.simulator.encounterDist        = 0.1;   % The distance for which two monomer are considered to have met 
            obj.simulator.showSimulation       = true; 
            obj.simulator.recordData           = false;
            obj.simulator.notifyByEmail        = false;
            obj.simulator.notifyCycleLength    = 30;    % number of simulation cycles after which an email is sent 
            obj.simulator.recipeFileName       = 'debugRecipe';
            obj.simulator.recipesFolder        = ''; 
        end
        
        function SetChainParams(obj)
            obj.chain = ChainParams;
        end
        
        function SetDomainParams(obj)
            obj.domain = DomainHandlerParams;
        end
        
        function SetDataRecorderParams(obj)
            obj.dataRecorder = SimulationDataRecorderParams;
        end
        
        function SetPlotHandlerParams(obj)
        end
        
    end
end