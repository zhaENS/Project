% scrAnalyzeTADData
% analyze the experimental data from luca et al

% class parameters
params.xlsFilePath = fullfile(pwd,'..','Data','Luca','MC_TAD_DE_E14d0_rep2+3_true.xlsx');
params.fillGapsBy   = 'sameValuesAsBoundary'; % when the data is  missing (i.e, the segment takes the position of few beads,
params.beadRangeToAnalyze.bead1 = 1:307; % if empty, analyze all,
params.beadRangeToAnalyze.bead2 = 1:307;

a = AnalyzeEncounterFrequencies;
a.Initialize(params)

