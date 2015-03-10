% scrRunReconstructionUsingRouse
% load data
load(fullfile(pwd,'..','..','PolymerChainDynamics','Code','ExperimentDataAnalysis','savedAnalysisTADDAndE.mat'));
br.bead1=1:307; br.bead2=1:307; [~, eMat2]= a.GetEncounterFrequencyMatrix(br,'average');

% load(fullfile(pwd,'..','..','PolymerChainDynamics','Documents\lab Meetings\09-01-15\simpleRouse_simpleRouse_TestDifferentPolymerStructuresForReconstruction_18_47_04-Jan-2015.mat'));
% load('D:\Ofir\Work\ENS\Polymer Chain Dynamics\Documents\lab Meetings\09-01-15\simpleRouse_simpleRouse_TestDifferentPolymerStructuresForReconstruction_18_47_04-Jan-2015.mat');
% eMat2=srf.beadEncounterHistogram.twoSides;
% get encounter frequency matrix

% initialize the constructor class
c=CalculateBeadDistancesByRouseModel;

% start the process
c.Initialize(eMat2);