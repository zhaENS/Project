% scrRunRouseSimulator
% function scrRunRouseSimulator
close all; 
close hidden
dbstop if error

% Initialize environment
% curPath = pwd;
% addpath(genpath(fullfile(curPath,'..','..','Utils')));
% xmlStr = fileread('SimulationFrameworkParams.xml');
% params = xml_parse(xmlStr);
simulatorParams = SimulationFrameworkParams('dt',0.01,'dimension',3,'numSteps',Inf);

% Initialize domains
domainForces    = ForceManagerParams('dt',simulatorParams.simulator.dt,'diffusionForce',true,...
                                     'diffusionConst',0.0001,'lennardJonesForce',false);
domainParams(1) = DomainHandlerParams('domainShape','sphere','domainWidth',10,'forceParams',domainForces);
% domainParams(2) = DomainHandlerParams('domainShape','cylinder','reflectionType','off',...
%                                       'diffusionForce',false,'domainWidth',1,'domainHeight',50);
                                  
% Initialize chains 
chainForces    = ForceManagerParams('dt',simulatorParams.simulator.dt,'springForce',true,...
                    'springConst',1,'bendingElasticityForce',false,'minParticleEqDistance',1,'bendingConst',1);
chainParams(1) = ChainParams('numBeads',128,'initializeInDomain',1,'forceParams',chainForces,'b',1);

chainParams(2) = ChainParams('numBeads',128,'initializeInDomain',1,'fixedBeadNum',[],'fixedBeadsPosition',[],'forceParams',chainForces);


chainParams(3) = ChainParams('numBeads',128,'initializeInDomain',1,'forceParams',chainForces);


chainParams(4) = ChainParams('numBeads',128,'initializeInDomain',1,'forceParams',chainForces);
                                     
% register parameters
simulatorParams.SetDomainParams(domainParams);
simulatorParams.SetChainParams(chainParams);

% SimulatorParams
% profile on 
r = RouseSimulatorFramework(simulatorParams);
r.Run
% profile viewer