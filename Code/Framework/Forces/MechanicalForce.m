classdef MechanicalForce<handle
    % This is a point force from forceCenter
    % estent is the radius of the force
    % drop is the force drop with distance
properties
    forceCenter
    extent
    
end

methods
    function obj = MechanicalForce
    end
    
    function newPos = Apply(obj,ParticlePosition)
    end
    
end

end
