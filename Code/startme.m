% Initialize environment
dbstop if error
curPath = pwd;
addpath(genpath(pwd));
addpath(genpath(fullfile(curPath,'..','..','Utils')));