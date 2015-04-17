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
params = SimulationFrameworkParams;

% SimulatorParams
profile on 
r = RouseSimulatorFramework(params);
r.Run
profile viewer