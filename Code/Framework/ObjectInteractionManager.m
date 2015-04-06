classdef ObjectInteractionManager<handle 
    % This class handles the interactions between objects in the domain 
    %TODO: finish 
    properties
        forceManager
    end
    
    events 
    end
    
    methods
        function obj = ObjectInteractionManager
            % assign force handler 
            obj.forceManager = ForceManager;
        end
    end
end