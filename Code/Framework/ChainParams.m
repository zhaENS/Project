classdef ChainParams<handle
    % this class holds the parameters for a chain Class 
    % TODO: consider creating individual class for each force
    
    properties
        dimension@double
        beta@double
        b@double
        diffusionConst@double
        dt@double
        numBeads@double
        connectedBeads@double
        bendingElasticityForce@logical
        lennardJonesForce@logical
        springForce@logical
        diffusionForce@logical
        minBeadDistance@double
        fixedBeadNum@double
        noiseDistribution@char
        noiseStd@double
        noiseMean@double
        springConst@double
        LJPotentialDepth@double
        LJPotentialWidth@double
        bendingConst@double
        stickyBeads@double
        allowSelfAffinity@logical
        forceParams
        
    end
    
    methods
        function obj = ChainParams                        
            obj.dimension              = 3;
            obj.beta                   = 2; % for rouse, place 2. 
            obj.b                      = sqrt(3);
            obj.dt                     = 1e-2;
            obj.diffusionConst         = 1;
            obj.numBeads               = 32;
            obj.connectedBeads         = [1 10; 1 15; 1 32];
            obj.bendingElasticityForce = false;
            obj.lennardJonesForce      = false;
            obj.springForce            = false;
            obj.diffusionForce         = false;
            obj.LJPotentialDepth       = 0.1;
            obj.LJPotentialWidth       = 0.1;
            obj.minBeadDistance        = 0;
            obj.fixedBeadNum           = [];
            obj.allowSelfAffinity      = false; % can sticky beads stick to other sticky beads on the same chain
            obj.stickyBeads            = [];    % beads that can stick to others, is used to stick to other chains  
            obj.noiseDistribution      = 'Gaussian';
            obj.noiseStd               = sqrt(2*obj.diffusionConst*obj.dt);
            obj.noiseMean              = 0;
            obj.springConst            = (obj.dimension*obj.diffusionConst./obj.b^2)*ones(obj.numBeads);
            obj.bendingConst           = 1; 
            
            % set spring constant for the connected beads
            for cIdx = 1:size(obj.connectedBeads,1)
                obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2))=obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2));
                obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1))=obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1));
            end
            
            %============================
            % Define force parameters 
            obj.forceParams                        = ForceManagerParams;
            
            obj.forceParams.bendingElasticityForce = obj.bendingElasticityForce;
            obj.forceParams.lennardJonesForce      = obj.lennardJonesForce;
            obj.forceParams.springForce            = obj.springForce;
            obj.forceParams.diffusionForce         = obj.diffusionForce;  
            
            obj.forceParams.LJPotentialDepth       = obj.LJPotentialDepth;
            obj.forceParams.LJPotentialWidth       = obj.LJPotentialWidth;
            obj.forceParams.bendingConst           = obj.bendingConst;
            obj.forceParams.dt                     = obj.dt;
            obj.forceParams.diffusionConst         = obj.diffusionConst;
            obj.forceParams.fixedParticleNum       = obj.fixedBeadNum;
            obj.forceParams.springConst            = obj.springConst;        
            obj.forceParams.minParticleDistance    = obj.minBeadDistance;
        end
        
    end
    
end