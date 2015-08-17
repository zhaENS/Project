classdef SimulationDataRecorderParams<handle
    % This class holds parameters for the simulation data recorder 
    properties
        saveType@char
        resultsFolder@char
        recipeFileName@char
        encounterDist@double
    end
    
    methods
        
        function obj = SimulationDataRecorderParams
            obj.SetDefaultParams
        end
        
        function SetDefaultParams(obj)
            obj.saveType       = 'external'; %save to external mat file. no data is save in the class. options [external/internal/none/all]
            %obj.resultsFolder  = 'D:\Ofir\Work\PolymerChainDynamicsResults';
           %  obj.resultsFolder =  'C:\projectsENS\PolymerChainDynamicsResults';
            obj.resultsFolder  = 'D:\Zha\Project\PolymerChainDynamicsResults';
            obj.recipeFileName = '';%will be copied by framework class
            obj.encounterDist  = []; %  The distance for which two monomer are considered to have met
        end
        
    end
end
