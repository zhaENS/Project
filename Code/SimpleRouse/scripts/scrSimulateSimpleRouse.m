% a script to run the simple rouse class. The exact procedure is described
% int the recipe file attached to each simulation. 

% load simulation parameters
rp = SimpleRouseParams;

% create the class
srf = SimpleRouseFramework; 

% start the simualtion
srf.Initialize(rp);

% save the class
d      = date;
c      = clock;
hour   = num2str(c(4));
minute = num2str(c(5));
resultsFolder = fullfile(pwd,'..','..','PolymerChainDynamicsResults');
saveName      = sprintf('%s',['simpleRouse_',rp.recipeFileName,'_',hour,'_', minute,'_',d,'.mat']);
% Remove bead distance and save
resultsStruct = struct;

% Remove bead distance and save 

if srf.params.saveResults
    srf.beadDistance = [];
    save(fullfile(resultsFolder,saveName),'srf')
end
