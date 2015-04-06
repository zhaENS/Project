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
%         lennardJonesForce@logical
        springForce@logical
%         diffusionForce@logical
        minBeadDistance@double
        fixedBeadNum@double
%         noiseDistribution@char
%         noiseStd@double
%         noiseMean@double
        springConst@double
%         LJPotentialDepth@double
%         LJPotentialWidth@double
        bendingConst@double
        stickyBeads@double
        allowSelfAffinity@logical
        forceParams
        
    end
    
    methods
        function obj = ChainParams                        
            obj.dimension              = 3;       % inherited from framework
            obj.beta                   = 2;       % for rouse, place 2. 
            obj.b                      = sqrt(3); % std of distance between beads
            obj.dt                     = 1e-2;    % inherited from framework
            obj.diffusionConst         = 1;
            obj.numBeads               = 32;
            obj.connectedBeads         = [1 10; 1 15; 1 32];
            obj.bendingElasticityForce = false;
            obj.springForce            = false;
            obj.minBeadDistance        = 0;
            obj.fixedBeadNum           = [];
            obj.allowSelfAffinity      = false; % can sticky beads stick to other sticky beads on the same chain
            obj.stickyBeads            = [];    % beads that can stick to others, is used to stick to other chains  
            obj.springConst            = (obj.dimension*obj.diffusionConst./obj.b^2)*ones(obj.numBeads);
            obj.bendingConst           = 1; 
            
            % set spring constant for the connected beads
            for cIdx = 1:size(obj.connectedBeads,1)
                obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2))=obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2));
                obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1))=obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1));
            end
            
            % Define force parameters 
            obj.forceParams  = ForceManagerParams;
            
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
            obj.forceParams.minParticleDistance    = obj.minBeadDistance;
            
        end
    end
    
end