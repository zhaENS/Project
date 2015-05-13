classdef ChainParams<handle
    % This class holds the parameters for a chain Class 
    
    properties
        dimension@double
        beta@double
        b@double
        diffusionConst@double
        dt@double
        numBeads@double
        connectedBeads@double
        bendingElasticityForce@logical
        springForce@logical
        minBeadDistance@double
        fixedBeadNum@double
        fixedBeadsPosition@double
        springConst@double
        bendingConst@double
        stickyBeads@double
        beadsOnBoundary@double
        allowSelfAffinity@logical
        forceParams        
    end
    
    methods
        
        function obj = ChainParams(varargin)
            % default 
            obj.dimension              = 3;         % inherited from framework
            obj.beta                   = 2;         % for rouse, place 2. 
            obj.b                      = 1*sqrt(3); % std of distance between beads
            obj.dt                     = 1e-2;      % inherited from framework
            obj.diffusionConst         = 1;
            obj.numBeads               = 32;
            obj.connectedBeads         = [];
            obj.bendingElasticityForce = false;
            obj.springForce            = true;           
            obj.fixedBeadNum           = [];    % beads which do not move 
            obj.fixedBeadsPosition     = [];    % position for the fixed beads (can be on the boundary)
            obj.beadsOnBoundary        = [];    % list of bead attached to the domain's boundary 
            obj.allowSelfAffinity      = false; % can sticky beads stick to other sticky beads on the same chain?
            obj.stickyBeads            = [];    % beads that can stick to others, is also used to stick to other chains                      
            obj.bendingConst           = 1;             
            
            if ~isempty(varargin)
                % parse input 
                 obj.ParseInputParams(varargin);             
            end
            
            obj.springConst = 1.0*(obj.dimension*obj.diffusionConst./obj.b^2)*ones(obj.numBeads);% defined as a matrix for all beads
            % set spring constant for the connected beads
            for cIdx = 1:size(obj.connectedBeads,1)
                obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2))=...
                    obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2));
                obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1))=...
                    obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1));
            end
             obj.minBeadDistance        = 0*ones(obj.numBeads);
             
            % Define force parameters 
            obj.forceParams  = ForceManagerParams;
            obj.SetForceParams;
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
            
            obj.forceParams.bendingElasticityForce = obj.bendingElasticityForce;
%             obj.forceParams.lennardJonesForce      = obj.lennardJonesForce;
            obj.forceParams.springForce            = obj.springForce;            
            obj.forceParams.bendingConst           = obj.bendingConst;
            obj.forceParams.dt                     = obj.dt;
            obj.forceParams.diffusionConst         = obj.diffusionConst;
            obj.forceParams.fixedParticleNum       = obj.fixedBeadNum;
            obj.forceParams.springConst            = obj.springConst;        
            obj.forceParams.minParticleEqDistance  = obj.minBeadDistance;
            
        end
    end
    
end