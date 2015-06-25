classdef DomainHandlerParams<handle
    % This class holds the parameters of the domain handler
    
    % TODO: delete parameters related to forces from the properties
    % TODO: add domain shapes: cube, polygon, and user defined/ import from meshlab/ blender
    % TODO: consider adding shape primatives for the geometry of the domain
    properties
        dimension@double
        domainShape@char % [sphere |cylinder |twoPlates |open] 
        domainWidth@double
        domainHeight@double
        domainLength@double
        LJPotentialWidth@double
        LJPotentialDepth@double
        lennardJonesForce@logical% lennard jones force on particles in the domain 
        diffusionForce@logical   % thermal (diffusion) force on particles in the domain 
        morseForce@logical
        morsePotentialDepth@double
        morsePotentialWidth@double
        morseForceType@char          % type of morse force [attractive|repulsive|full]
        minParticleEqDistance@double % minimum equilibrium distance 
        diffusionConst@double    % diffusion constant 
        reflectionType@char      % current only option: [preserveEnergy | off]
        reflectionDirection@char % either [in |out]
        domainCenter@double   % position in space for the center default [0 0 0] 3d
        showDomain@logical    % flag to show or hide the domain 
        dt@double             % simulation step time 
        permeability@double   % how permiable is the domain surface range:[0-1] % to be used in the future 
        maxReflectionsPerParticle@double % the maximal number of times a particle can be reflected in a single reflection loop
        parentAxes % obsolete
        forceParams = ForceManagerParams;
    end
    
    methods
         
        function obj = DomainHandlerParams(varargin)
            
            % Default parameters             
            obj.domainShape                = 'sphere'; %[sphere | cylinder |twoplates| open]
            obj.dimension                  = 3;
            obj.domainWidth                = 15;
            obj.domainHeight               = 15;
            obj.domainCenter               = [0 0 0];
            obj.showDomain                 = true; % obsolete- move to graphics
            obj.reflectionType             = 'preserveEnergy'; % reflection type 
            obj.reflectionDirection        = 'in'; % direction of reflection
            obj.maxReflectionsPerParticle  = 100;
            obj.permeability               = 0; % how permeable is the surface (saved for future use)  
                        
            % set input parameters            
            obj.ParseInputParams(varargin);
            
        end
        
        function ParseInputParams(obj,varargin)
            % parse name-value pair input parameters
            v = varargin{:};
            if mod(numel(varargin{:}),2)~=0
                error('name-value pair argument must be inserted')
            end
            
            for vIdx = 1:(numel(v)/2)
                obj.(v{2*vIdx-1})= v{2*vIdx};
            end
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
