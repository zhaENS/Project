% scrRunRouseSimulator
% function scrRunRouseSimulator
close all; 
close hidden
% clear hidden
% clear all
% clear classes
dbstop if error
profile on 

% Initialize environment
curPath = pwd;
addpath(genpath(fullfile(curPath,'..','..','Utils')));
xmlStr = fileread('SimulationFrameworkParams.xml');
params = xml_parse(xmlStr);

% SimulatorParams
r = RouseSimulatorFramework(params);
r.Run
