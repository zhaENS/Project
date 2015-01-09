% scrTestReconstructionByRouse
% for the experimental data
c=CalculateBeadDistancesByRouseModel;
c.smoothingSpan = 2;
c.beadRange      = struct('bead1',1:307,...
                        'bead2',1:307);
                    c.Initialize;
% for the test structure
% load the data 
load(fullfile(pwd,'..','..','PolymerChainDynamics','Documents\lab Meetings\09-01-15\simpleRouse_simpleRouse_TestDifferentPolymerStructuresForReconstruction_18_47_04-Jan-2015.mat'));
eMat = srf.beadEncounterProbability.oneSide;

c2=CalculateBeadDistancesByRouseModel;
c2.smoothingSpan = 5;
c2.beadRange      = struct('bead1',1:128,...
                        'bead2',1:128);
c2.Initialize(eMat);
                    
                    