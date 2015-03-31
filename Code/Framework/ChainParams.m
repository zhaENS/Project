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
        minBeadDist@double
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
        
    end
    
    methods
        function obj = ChainParams
            obj.SetDefaultParams;
        end
                
        function SetDefaultParams(obj)
            
            obj.dimension              = 3;
            obj.beta                   = 2; % for rouse, place 2. 
            obj.b                      = sqrt(3);
            obj.dt                     = 1e-2;
            obj.diffusionConst         = 1;
            obj.numBeads               = 64;
            obj.connectedBeads         = [];
            obj.bendingElasticityForce = false;
            obj.lennardJonesForce      = false;
            obj.springForce            = true;
            obj.diffusionForce         = true;
            obj.minBeadDist            = 0;
            obj.fixedBeadNum           = [];
            obj.allowSelfAffinity      = false; % can sticky beads stick to other sticky beads on the same chain
            obj.stickyBeads            = [];    % beads that can stick to others, is used to stick to other chains  
            obj.noiseDistribution      = 'Gaussian';
            obj.noiseStd               = sqrt(2*obj.diffusionConst*obj.dt);
            obj.noiseMean              = 0;
            obj.springConst            = (obj.dimension*obj.diffusionConst./obj.b^2)*ones(obj.numBeads);
            % set spring constant for the connected beads
            for cIdx = 1:size(obj.connectedBeads,1)
                obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2))=obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2));
                obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1))=obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1));
            end
            
            obj.LJPotentialDepth = 10^(-2);
            obj.LJPotentialWidth = 10^(-2);
            obj.bendingConst     = 90;
        end
        
    end
    
end