classdef ForceManager<handle
    % This class collects different forces acting upon a body/particle
    % system.
    % Forces are considered to be addative, hence the order of their
    % operations is not important
    % forceManager is a child of a physical object like a chain or a
    % particle system
    % connectivityMap should be logical
    % TODO: add a listener to object parameter change
    % TODO: add an option to recreate the noise values from previously used
    %       seed or simulation parameters
    properties
        %  Names of forces acting on the particle system
        params
        springForce            = false;
        lennardJonesForce      = false;
        diffusionForce         = false;
        bendingElasticityForce = false;
        edges % matrices representing the edges between connected particles
        particleDistance % pairwise distance between particles
        
    end
    
    methods
        
        function obj = ForceManager
            
        end
        
        function newParticlePosition = Apply(obj,particlePosition,connectivityMap,dt,...
                diffusionConst,springConst,LJPotentialWidth,LJPotentialDepth,bendingConst,minParticleDist,fixedParticleNum)
            % Apply chosen forces on the N vertices/particles of an object
            % object must contain the position of its vertices in any
            % dimension
            % particlePosition should be an NxDim array of coordinates
            % particledist is the pairwise distance between the particles
            % NxN matrix
            % connectivityMap is a binary graph representing which particle
            % is connected to which other particle
            
            % Calculate the edges vectors in all dimension
            obj.edges   = GetEdgesVectors_mex(particlePosition,connectivityMap);
            
            % get pair-wise distance between particles
            particleDistances = obj.GetParticleDistance(particlePosition);
            obj.particleDistance = particleDistances;% save property
            
            % obtain the edges length
            edgeLength = obj.GetEdgesLength(particleDistances,connectivityMap);
            
            % Apply forces
            springForces =  obj.GetSpringForce(edgeLength,springConst,connectivityMap,minParticleDist,fixedParticleNum);
            
            % Lenard-jones force
            ljForces      = obj.GetLenardJonesForce(particlePosition,particleDistances,LJPotentialWidth,LJPotentialDepth);
            
            % thermal (diffusion) force
            diffusionForces = obj.GetDiffusionForce(particlePosition,diffusionConst,dt);
            
            % bending forces
            bendingForces = obj.GetBendingElasticityForce(obj.edges,connectivityMap,bendingConst);
            % the effect of appliying forces is considered addative
            newParticlePosition = -springForces*dt*particlePosition+...
                ljForces*dt+...
                bendingForces*dt+...
                diffusionForces+...
                particlePosition;
        end
        
        function force  = GetSpringForce(obj,particleDist,springConst,connectivityMap,minParticleDist,fixedParticleNum)
            % Calculate the spring force between N particles in any dimension M.
            % particleDist    - NxN matrix of pairwise particle distances
            % springConst     - NxN double matrix of spring constants
            % connectivityMap - NxN binary matrix which defines the connectivity between particles
            % minParticleDist - minimal distance between particles
            % fixedParticleNum - particles in the system which do not move
            force = zeros(size(connectivityMap));
            if obj.springForce
                connectivityMap           = (connectivityMap);
                L                         = (1-minParticleDist./particleDist).*connectivityMap;
                L(~connectivityMap)       = 0;
                sumForces                 = sum(L,2);
                force                     = springConst.*(diag(sumForces)-L); % set the maindiagonal
                force(fixedParticleNum,:) = 0;% zero out forces for fixed particles
                force(:,fixedParticleNum) = 0;% zero out forces for fixed particles
            end
        end
        
        function force  = GetDiffusionForce(obj,particlePosition,diffusionConst,dt)
            % get thermal (diffusion) force
            force = zeros(size(particlePosition));
            if obj.diffusionForce
                force = randn(size(particlePosition))*sqrt(2*diffusionConst*dt);
            end
            
        end
        
        function force  = GetLenardJonesForce(obj,particlePosition,particleDist,LJPotentialWidth,LJPotentialDepth)% move to mex
            % calculate Lenard jones force between particles
            force = zeros(size(particlePosition));
            if obj.lennardJonesForce
                force = LennardJones_mex(particlePosition,particleDist,LJPotentialWidth,LJPotentialDepth);
            end
        end
        
        function force  = GetBendingElasticityForce(obj,edgeMat,connectivityMat,bendingConst)% needs clean up
            
            force = zeros(size(connectivityMat,1),size(edgeMat,3));
            if obj.bendingElasticityForce
                force = BendingElasticityForce_mex(edgeMat(:,:,1),edgeMat(:,:,2), edgeMat(:,:,3), connectivityMat,bendingConst);
            end
            
        end
    end
    
    methods (Static)
        
        function edgesVec   = GetEdgesVectors(particlePosition, connectivityMap)
            % Calculate the edges vectors between connected particles
            edgesVec = zeros(size(particlePosition,1),size(particlePosition,1), size(particlePosition,2));
            for dIdx=1:size(particlePosition,2)
                edgesVec(:,:,dIdx) = bsxfun(@times,particlePosition(:,dIdx),connectivityMap)-...
                    bsxfun(@times,particlePosition(:,dIdx)',connectivityMap);
            end
        end
        
        function edgeLength = GetEdgesLength(edgeVec,connectivityMap)
            % Calculate the lengths of edges
            edgeLength = edgeVec.*connectivityMap;
        end
        
        function particleDistance = GetParticleDistance(particlePosition)
            % Get pair-wise particel distance using the mex function
            % pdist2mex
            particleDistance = pdist2mex(particlePosition',...
                particlePosition','euc',[],[],[]);
        end
        
    end
    
end
