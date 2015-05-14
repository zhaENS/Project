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


% Initialize domains
domainParams(1) = DomainHandlerParams('domainShape','sphere','domainWidth',10,'diffusionForce',true);
domainParams(2) = DomainHandlerParams('domainShape','cylinder','reflectionType','off',...
                                      'diffusionForce',false,'domainWidth',1,'domainHeight',50);
                                  
% Initialize chains 
chainParams(1) = ChainParams('numBeads',10,'initializeInDomain',1);
% chainParams(1).fixedBeadNum = 5;
% chainParams(1).fixedBeadsPosition = randn(1,3);

chainParams(2) = ChainParams('numBeads',100,'initializeInDomain',1,'fixedBeadNum',[],'fixedBeadsPosition',[]);
% chainParams(2).fixedBeadNum = [5 10 15 20 25 30];
% chainParams(2).fixedBeadsPosition = randn(6,3);

chainParams(3) = ChainParams('numBeads',64,'initializeInDomain',1);
% chainParams(3).fixedBeadNum = [1 4 12 27 32 56];
% chainParams(3).fixedBeadsPosition = randn(6,3);

chainParams(4) = ChainParams('numBeads',32','initializeInDomain',1);
% chainParams(4).fixedBeadNum = [5 10 15 20 25 30];
% chainParams(4).fixedBeadsPosition = randn(6,3);
                                     
params = SimulationFrameworkParams(chainParams,domainParams);

% SimulatorParams
% profile on 
r = RouseSimulatorFramework(params);
r.Run
% profile viewer