classdef DomainHandlerParams<handle
    % this class holds the parameters of the domain handler
    properties
        dimension@double
        domainShape@char % sphere/cylinder/twoPlates/none TODO: add cube, add user defined/ import from meshlab
        domainWidth@double
        domainHeight@double
        LJPotentialWidth@double
        LJPotentialDepth@double
        lennardJonesForce@logical% lennard jones force on particles in the domain 
        diffusionForce@logical % thermal (diffusion) force on particles in the domain 
%         bendingElasticityForce@logical
        diffusionConst@double % diffusion constant 
%         springForce@logical 
        reflectionType@char % current only option: preserveEnergy
        domainCenter@double % position in space for the center default [0 0 0] 3d
        showDomain@logical  % flag to show or hide the domain 
        dt@double           % simulation step time 
        maxReflectionsPerParticle@double % the maximal number of times a particle can be reflected in a single reflection loop
        parentAxes
        forceParams
    end
    
    methods
        function obj = DomainHandlerParams
            obj.domainShape            = 'sphere';
            obj.domainWidth            = 20;
            obj.domainHeight           = 10;
            obj.domainCenter           = [0 0 0];
            obj.showDomain             = true;
            
            % forces of the domain 
            obj.LJPotentialWidth       = 0.5;
            obj.LJPotentialDepth       = 0.1;
            obj.diffusionConst         = 1;
            obj.lennardJonesForce      = false;
            obj.diffusionForce         = true;
            
            obj.dt                     = 0;%  assigned by Framework
            obj.reflectionType         = 'preserveEnergy';
            obj.maxReflectionsPerParticle  = 30;
                        
            obj.forceParams                = ForceManagerParams;
            obj.SetForceParams;
        end
        
        function SetForceParams(obj)
            % forces inherit class params
            obj.forceParams.dt                     = obj.dt;
            obj.forceParams.lennardJonesForce      = obj.lennardJonesForce;
            obj.forceParams.diffusionConst         = obj.diffusionConst;
            obj.forceParams.diffusionForce         = obj.diffusionForce;
%             obj.forceParams.bendingElasticityForce = obj.bendingElasticityForce;
%             obj.forceParams.springForce            = obj.springForce;
            obj.forceParams.LJPotentialDepth       = obj.LJPotentialDepth;
            obj.forceParams.LJPotentialWidth       = obj.LJPotentialWidth;
        end
    end
end
