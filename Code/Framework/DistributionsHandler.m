classdef DistributionsHandler<handle
    % calculate the distributions related to quantities of the rouse
    % polymer model 
    properties
        handles
        numChains
        distributions
        frequencies
        positions
        times
        step
        simulationRound
    end
    
    events
    end
    
    methods
        function obj = DistributionsHandler(rouseObj)
            % construct the DistributionHandler with Rouse chain objects;
            obj.step                  = 0;  
            obj.simulationRound       = 1;
            obj.numChains             = numel(rouseObj);
            obj.handles.classes.rouse = rouseObj;
            for cIdx = 1:obj.numChains
              obj.frequencies.beadDist.chain{cIdx}                    = zeros(obj.handles.classes.rouse(cIdx).params.numBeads);
              obj.distributions.beadDist.chain(cIdx).mean             = zeros(obj.handles.classes.rouse(cIdx).params.numBeads);
              obj.distributions.beadDist.chain(cIdx).adjacentBeadsRMS = 0;
            end
                     
            % For all chains 
%             obj.frequencies.beadDist.total = zeros(obj.handles.classes.rouse(cIdx).params.numBeads.value);
        end 
        
        function Add(obj)
            % add data to the stored distributions
            obj.step = obj.step+1;% each call to Add raises the step count
            obj.AddFrequencies
            obj.AddPositions
%             obj.CalcDistributions
        end
        
        function AddFrequencies(obj)
            for cIdx = 1:obj.numChains
             obj.frequencies.beadDist.chain{cIdx}(:,:,obj.step) = obj.handles.classes.rouse(cIdx).beadsDist;
            end
        end        
        
        function AddPositions(obj)
            % Record the position of the rouse chains 
            for cIdx = 1:obj.numChains
                obj.positions.chains{cIdx}(obj.step) = obj.handles.classes.rouse(cIdx).positions.beads.cur;
            end            
        end
        
        function CalcDistributions(obj)%TODO: check
                for cIdx = 1:obj.numChains
                  obj.distributions.beadDist.chain(cIdx).mean  = mean(obj.frequencies.beadDist.chain{cIdx},3);  
                  obj.distributions.beadDist.chain(cIdx).std   = std(obj.frequencies.beadDist.chain{cIdx},0,3);
                   
                  s = sqrt(mean(obj.frequencies.beadDist.chain{cIdx},3).^2);
                  s = s(triu(obj.handles.classes.rouse(cIdx).connectionMap.indices.in.map));

                  obj.distributions.beadDist.chain(cIdx).adjacentBeadsRMS = sqrt(sum(s)/(numel(s)*obj.step));
                      
                  obj.distributions.beadDist.endToEnd  = mean(obj.frequencies.beadDist.chain{cIdx}(1,end,:).^2);

                end
        end     
        
        function NewSimulation(obj)
            
        end
    end
end
