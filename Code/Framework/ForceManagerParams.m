classdef ForceManagerParams<handle
    % parameter class for the forceManger class
    properties
        springForce            = false;
        lennardJonesForce      = false;
        diffusionForce         = false;
        bendingElasticityForce = false;
        morseForce             = false;
        % Parameters for the forces
        bendingConst           = 0;
        springConst            = 0;
        diffusionConst         = 0;
        LJPotentialWidth       = 0;
        LJPotentialDepth       = 0;
        minParticleEqDistance  = 0;
        dt                     = 0; % the time to activate the force
        fixedParticleNum       = [];
        morsePotentialDepth    = 0;
        morsePotentialWidth    = 0;
        morseForceType         = '';
        edges % matrices representing the edges between connected particles
        particleDistance % pairwise distance between particles
    end 
    
    methods 
        
        function obj = ForceManagerParams(varargin)
                % parse name-value pair input parameters
                if ~isempty(varargin)
                    v = varargin;
                    if mod(numel(varargin{:}),2)~=0
                        error('name-value pair argument must be inserted')
                    end
                    
                    for vIdx = 1:(numel(v)/2)
                        obj.(v{2*vIdx-1})= v{2*vIdx};
                    end
                end
        end
    end
end