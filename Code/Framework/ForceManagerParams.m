classdef ForceManagerParams<handle
    % parameter class for the forceManger class
    %TODO: make a class for each force
    properties
        % force type
        springForce@logical
        lennardJonesForce@logical
        diffusionForce@logical
        bendingElasticityForce@logical
        morseForce@logical
        mechanicalForce@logical % mechanical force acting on particles (push/pull/etc)
        
        % Parameters for the forces
        bendingConst@double
        openningAngle@double
        springConst@double
        diffusionConst@double
        LJPotentialWidth@double
        LJPotentialDepth@double
        minParticleEqDistance@double % can be a matrix representing pairwise min particle distance
        dt@double % the time to activate the force [obsolete]
        fixedParticleNum@double
        morsePotentialDepth@double
        morsePotentialWidth@double
        morseForceType@char
        mechanicalForceCenter@double % point source center
        mechanicalForceDirection@char
        mechanicalForceMagnitude@double
        edges@double %matrices representing the edges between connected particles
        particleDistance@double % pairwise distance between particles (unused)
    end 
    
    methods 
        
        function obj = ForceManagerParams(varargin)
                % parse name-value pair input parameters
                if ~isempty(varargin)
                    v = varargin;
                    if mod(numel(varargin),2)~=0
                        error('name-value pair argument must be inserted')
                    end
                    
                    for vIdx = 1:(numel(v)/2)
                        obj.(v{2*vIdx-1})= v{2*vIdx};
                    end
                end
        end
    end
end