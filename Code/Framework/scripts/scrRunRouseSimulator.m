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

% initialize chain parameters
chainParams(1) = ChainParams('numBeads',10);
chainParams(2) = ChainParams('numBeads',100);
chainParams(3) = ChainParams('numBeads',64');

params = SimulationFrameworkParams(chainParams);

% SimulatorParams
% profile on 
r = RouseSimulatorFramework(params);
r.Run
% profile viewer