classdef DomainHandlerParams<handle
    % this class holds the parameters of the domain handler
    properties
        dimension@double
        domainShape@char % sphere/cylinder/twoPlates/open TODO: add cube, add user defined/ import from meshlab
        domainWidth@double
        domainHeight@double
        LJPotentialWidth@double
        LJPotentialDepth@double
        lennardJonesForce@logical% lennard jones force on particles in the domain 
        diffusionForce@logical   % thermal (diffusion) force on particles in the domain 
        morseForce@logical
        morsePotentialDepth@double
        morsePotentialWidth@double
        morseForceType@char          % type of morse force [attractive|repulsive|full]
        minParticleEqDistance@double % minimum equilibrium distance 
        diffusionConst@double % diffusion constant 
        reflectionType@char   % current only option: preserveEnergy
        domainCenter@double   % position in space for the center default [0 0 0] 3d
        showDomain@logical    % flag to show or hide the domain 
        dt@double             % simulation step time 
        maxReflectionsPerParticle@double % the maximal number of times a particle can be reflected in a single reflection loop
        parentAxes
        forceParams
    end
    
    methods
         
        function obj = DomainHandlerParams
            obj.domainShape            = 'sphere'; %[sphere | cylinder |twoplates| open
            obj.domainWidth            = 5;
            obj.domainHeight           = 5;
            obj.domainCenter           = [0 0 0];
            obj.showDomain             = true; % move to graphics
            
            % forces of the domain 
            obj.LJPotentialWidth       = 0.2;
            obj.LJPotentialDepth       = 0.2;
            obj.diffusionConst         = 0;   % assigned by framework
            obj.lennardJonesForce      = false;
            obj.diffusionForce         = true;   
            obj.morseForce             = false; 
            obj.morsePotentialDepth    = .01;
            obj.morsePotentialWidth    = .01;
            obj.morseForceType         = 'attractive'; % type of morse force [attractive|repulsive|full]
            obj.minParticleEqDistance  = 0.2; % min particle equilibrium distance
            obj.dt                     = 0;   % assigned by Framework
            obj.reflectionType         = 'preserveEnergy';
            obj.maxReflectionsPerParticle  = 100;
                        
            obj.forceParams                = ForceManagerParams;
            obj.SetForceParams;
        end
        
        function SetForceParams(obj)
            % forces inherit class params
            obj.forceParams.dt                     = obj.dt;
            obj.forceParams.lennardJonesForce      = obj.lennardJonesForce;     
            obj.forceParams.morseForce             = obj.morseForce;
            obj.forceParams.LJPotentialDepth       = obj.LJPotentialDepth;
            obj.forceParams.LJPotentialWidth       = obj.LJPotentialWidth;
            obj.forceParams.morsePotentialDepth    = obj.morsePotentialDepth;
            obj.forceParams.morsePotentialWidth    = obj.morsePotentialWidth;
            obj.forceParams.morseForceType         = obj.morseForceType;
            obj.forceParams.minParticleEqDistance  = obj.minParticleEqDistance;
            obj.forceParams.diffusionConst         = obj.diffusionConst;
            obj.forceParams.diffusionForce         = obj.diffusionForce;
        end
    end
end
