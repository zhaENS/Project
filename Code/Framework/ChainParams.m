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
        minBeadDist@double
        fixedBeadNum@double
        noiseDistribution@char
        noiseStd@double
        noiseMean@double
        springConst@double
        LJPotentialDepth@double
        LJPotentialWidth@double
        bendingConst@double
        
    end
    
    methods
        function obj = ChainParams
            obj.SetDefaultParams;
        end
                
        function SetDefaultParams(obj)
            
            obj.dimension              = 3;
            obj.beta                   = 2; % for rouse, place 2. 
            obj.b                      = 0.1;
            obj.dt                     = 0.2;
            obj.diffusionConst         = 1;
            obj.numBeads               = 64;
            obj.connectedBeads         = [];
            obj.bendingElasticityForce = false;
            obj.lennardJonesForce      = false;
            obj.springForce            = true;
            obj.minBeadDist            = 0;
            obj.fixedBeadNum           =[];
            obj.noiseDistribution      = 'Gaussian';
            obj.noiseStd               = sqrt(2*obj.diffusionConst);
            obj.noiseMean              = 0;
            
            % set spring constant for the connected beads
            for cIdx = 1:size(obj.connectedBeads,1)
                obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2))=obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2));
                obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1))=obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1));
            end
            
            obj.LJPotentialDepth = 10^(-15);
            obj.LJPotentialWidth = 10^(-5);
            obj.bendingConst     = 90;
        end
        
    end
    
end