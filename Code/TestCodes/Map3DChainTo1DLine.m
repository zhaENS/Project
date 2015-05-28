function Map3DChainTo1DLine

close all 

% Define a simulationFramework params 
simulationParams = SimulationFrameworkParams('dt', 0.001,'numSteps',100,'dimension',3,...
    'objectInteraction',false,'showSimulation',true);

% define a domain 
domainForces     = ForceManagerParams('diffusionForce',true,'diffusionConst',0.01,'dt',simulationParams.simulator.dt);
openDomainParams = DomainHandlerParams('domainShape','open','forceParams',domainForces);

% define a chain
chainForces = ForceManagerParams('dt',simulationParams.simulator.dt,'springForce',true,...
    'springConst',1,'minParticleEqDistance',0);
chainParams = ChainParams('b',1,'numBeads', 50,'forceParams',chainForces);

% set the parameters in the simulation framework
simulationParams.SetChainParams(chainParams);
simulationParams.SetDomainParams(openDomainParams);

% initialize the framework 
simFW = RouseSimulatorFramework(simulationParams);
%------
simFW.PreSimulationBatchActions;
simFW.PreRunActions;
%------  

for sIdx = 1:simulationParams.simulator.numSteps
    simFW.Step
    
    % get the chain position 
    chainPos = simFW.objectManager.GetPosition(1);
    chainPos = chainPos{1};    
    % find the chain length
    chainLength = sum(sqrt(sum(diff(chainPos).^2,2)));
    
    
end




end
