classdef ForceManagerParams<handle
    % parameter class for the forceManger class
    properties
        springForce            = false;
        lennardJonesForce      = false;
        diffusionForce         = false;
        bendingElasticityForce = false;
        
        % Parameters for the forces
        bendingConst           = 0;
        springConst            = 0;
        diffusionConst         = 0;
        LJPotentialWidth       = 0;
        LJPotentialDepth       = 0;
        minParticleDistance    = 0;
        dt                     = 0; % the time to activate the force
        fixedParticleNum       = [];
                
        edges % matrices representing the edges between connected particles
        particleDistance % pairwise distance between particles
    end    
end