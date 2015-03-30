classdef Spring<handle
    % This class represents an harmonic spring object in polymer chain
    % a spring has head and tail points 
    properties
        springConst
        length
        minLength
    end
    
    events
    end
    
    methods
        
        function obj = Spring(parentObj)
            addlistener(parentObj,'getForce',@obj.CalculateForce)
        end
        
        function CalculateForce(obj,src,varargin)
            % Calculate the spring force between N particles in any dimension M.
            % particleDist    - NxN matrix of pairwise particle distances
            % springConst     - NxN double matrix of spring constants
            % connectivityMap - NxN binary matrix which defines the connectivity between particles
            % minParticleDist - minimal distance between particles
            % fixedParticleNum - particles in the system which do not move
            
            connectivityMap          = double(src.connectionMap.map);
            connectivityMap(isnan(connectivityMap)) = 0;
            % force              = springConst.*particleDist.*connectivityMap;
            particleDist = src.beadsDist+diag(Inf*ones(1,size(src.beadsDist,1)));
            L                         = (1-src.params.minBeadDist./particleDist).*connectivityMap;
            sumForces                 = sum(L,2);
            src.forces.spring         = src.params.springConst.*(diag(sumForces)-L); % set the maindiagonal
            src.forces.spring(src.params.fixedBeadNum,:) = 0;% zero out forces for fixed particles
            src.forces.spring(:,src.params.fixedBeadNum) = 0;% zero out forces for fixed particles
        end
        
    end
        
end
