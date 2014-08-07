classdef AnalyzeBeadConnectionSimulations<handle 
    % analyze the results of simulation of the Rouse polymer with N beads
    %  in each simulation bead j is connected to bead k, j=1:29, k=k
    %  j+2:end. This is a dedicated class for analyzing the results of the
    %  specific simulation.
    % Reminder: 
    % each batch contains simulation where bead j is connected to bead k,
    % where k= j+2:numBeads
    % there are 29 batches since connecting bead 31 to 33 is
    % symmetric as to connect bead 29 to 31. 
    properties
        resultsFolderPath
        numFiles
        fileNames
        filePath = {};
        simulationResults 
    end
    
    methods 
        function obj = AnalyzeBeadConnectionSimulations(resultsFolder)
            obj.resultsFolderPath = resultsFolder;
            obj.GetFileNames
            obj.Analyze
        end
        
        function GetFileNames(obj)
            % save file names and absolute path 
            d = dir([obj.resultsFolderPath,'\*.mat']);
            obj.numFiles = numel(d);
            for dIdx = 1:obj.numFiles
                obj.filePath{dIdx}  = fullfile(obj.resultsFolderPath,d(dIdx).name);
                obj.fileNames{dIdx} = d(dIdx).name;                
            end
        end
        
        function Analyze(obj)
            % calculate the eigenvalues and eigenvectors of the distance
            % matrices
            for dIdx=1:obj.numFiles
              [inds] = regexp(obj.fileNames{dIdx},'_');
              simulationBatchNumber =  str2double(obj.fileNames{dIdx}(inds(1)+1:inds(2)-1));
              load(obj.filePath{dIdx})
              connectedMonomers = results.params.rouseParams.connectMonomers;
              [v,e] = eig(results.beadDistMat,'nobalance');
              % TODO: complete the missing eigenvalues from previous
              % simulations 
              obj.simulationResults.batch(simulationBatchNumber).eigValues(:,connectedMonomers(2)-connectedMonomers(1)-1) = diag(e);
              obj.simulationResults.batch(simulationBatchNumber).eigVectors{connectedMonomers(2)-connectedMonomers(1)-1} = v;
              obj.simulationResults.batch(simulationBatchNumber).connectedMonomers(connectedMonomers(2)-connectedMonomers(1)-1,:) = connectedMonomers;
              % arrange the eig according to the simulated distance between
              % connected bead. i,e in beadDist(2) the eigenvalues from all
              % batches (in columns) for which the distance is 2 between
              % connected monomeners. distance is in the sense of index
              % difference
              obj.simulationResults.beadDist(connectedMonomers(2)-connectedMonomers(1)).eig(:,simulationBatchNumber) =diag(e);
              obj.simulationResults.beadEncounterHist(:,:,dIdx) = results.beadEncounterHist;
              % save the mean distance matrix 
              obj.simulationResults.beadDistMatrix(:,:,dIdx) = results.beadDistMat;
            end

        end
        
        function CalculateEigVariation(obj)% Unfinished
            % calculate the variation in the eigen values in bathes as
            % throughout the simulation of this batch. 
            for dIdx=1:obj.numFiles
              [inds] = regexp(obj.fileNames{dIdx},'_');
              simulationBatchNumber =  str2double(obj.fileNames{dIdx}(inds(1)+1:inds(2)-1));            
              obj.simulationResults.batch(simulationBatchNumber).eigValues(:,connectedMonomers(2)-connectedMonomers(1)-1) = diag(e);
              obj.simulationResults.batch(simulationBatchNumber).eigVectors{connectedMonomers(2)-connectedMonomers(1)-1} = v;
              obj.simulationResults.batch(simulationBatchNumber).connectedMonomers(connectedMonomers(2)-connectedMonomers(1)-1,:) = connectedMonomers;
              
            end
        end
    end
end