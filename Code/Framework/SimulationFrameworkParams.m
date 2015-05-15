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
        
        function obj = SimulationFrameworkParams()
            % class constructor 
            obj.SetSimulatorParams;
%             obj.SetChainParams(chainParams);
%             obj.SetDomainParams(domainParams);
            obj.SetDataRecorderParams;
            obj.SetPlotHandlerParams;
        end
        
        function SetSimulatorParams(obj)
            % Simulator framework params
            obj.simulator.dimension            = 3;
            obj.simulator.runSimulation        = false; % a flag indicating whether to allow the simulation to run at initiation 
            obj.simulator.numSimulationBatches = 1;     % number of simulation batches
            obj.simulator.numSimulations       = 1;     % number of simulations in each batch
            obj.simulator.numSteps             = 1000;   % for inf place Inf
            obj.simulator.dt                   = 1e-2;  % time step 
            obj.simulator.numChains            = 1;     % TODO: should be determined after chain initialization 
            obj.simulator.encounterDist        = 0.1;   % TODO: move out- The distance for which two monomer are considered to have met 
            obj.simulator.showSimulation       = true; 
            obj.simulator.recordData           = false;
            obj.simulator.notifyByEmail        = false;
            obj.simulator.notifyCycleLength    = 32;    % number of simulation cycles after which an email is sent 
            obj.simulator.recipeFileName       = 'debugRecipe';
            obj.simulator.recipesFolder        = ''; 
            
            % Control domain forces
            obj.simulator.diffusionConst       = 0.1; %TODO: should be determined for each domain seperately
%             obj.simulator.LJPotentialWidth     = 0.3;
%             obj.simulator.LJPotentialDepth     = 0.3;
        end
        
        function SetChainParams(obj,chainParams)%TODO: change force params accordingly
            % Assign parameters for each chain 
            % chain params should be a vector of ChainParams structures 
            
            if exist('chainParams','var')
                % change te simulation numChains parmeter
                obj.simulator.numChains = numel(chainParams);                
                obj.chain               = chainParams;                
            else
                
             for cIdx = 1:obj.simulator.numChains
                obj.chain(cIdx) = ChainParams;
             end
             
            end
            
            for cIdx = 1:obj.simulator.numChains
                % inherit the framework parameters 
                obj.chain(cIdx).forceParams.diffusionConst = obj.simulator.diffusionConst;
                obj.chain(cIdx).dimension                  = obj.simulator.dimension;            
                obj.chain(cIdx).forceParams.dt             = obj.simulator.dt;
%                 obj.chain(cIdx).SetForceParams;
            end
        end
        
        function SetDomainParams(obj,domainParams) 
            
            obj.domain = domainParams;
            
            % Inherit framework parameters
            for dIdx = 1:numel(obj.domain)
                obj.domain(dIdx).dimension        = obj.simulator.dimension;
                obj.domain(dIdx).forceParams.diffusionConst   = obj.simulator.diffusionConst;
                obj.domain(dIdx).showDomain       = obj.simulator.showSimulation;
                obj.domain(dIdx).forceParams.dt   = obj.simulator.dt;
                %             obj.domain.LJPotentialWidth = obj.simulator.LJPotentialWidth;
                %             obj.domain.LJPotentialDepth = obj.simulator.LJPotentialDepth;
%                 obj.domain(dIdx).SetForceParams;% update the force parameters of the domain 
            end
           
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