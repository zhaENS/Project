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

% Initialize chains 
chainParams(1) = ChainParams('numBeads',10);
chainParams(1).fixedBeadNum = 5;
chainParams(1).fixedBeadsPosition = randn(1,3);

chainParams(2) = ChainParams('numBeads',100);
chainParams(2).fixedBeadNum = [5 10 15 20 25 30];
chainParams(2).fixedBeadsPosition = randn(6,3);

chainParams(3) = ChainParams('numBeads',64');
chainParams(3).fixedBeadNum = [1 4 12 27 32 56];
chainParams(3).fixedBeadsPosition = randn(6,3);

chainParams(4) = ChainParams('numBeads',32');
chainParams(4).fixedBeadNum = [5 10 15 20 25 30];
chainParams(4).fixedBeadsPosition = randn(6,3);

params = SimulationFrameworkParams(chainParams);

% SimulatorParams
% profile on 
r = RouseSimulatorFramework(params);
r.Run
% profile viewer