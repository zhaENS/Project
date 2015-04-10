classdef ForceManager<handle
    % This class collects different forces acting upon a body/particle
    % system.
    % Forces are considered to be addative, hence the order of their
    % operations is not important
    % forceManager is a child of a physical object like a chain or a
    % particle system
    % connectivityMap should be logical
    % TODO: add an option to recreate the noise values from previously used
    %       seed or simulation parameters
    %TODO: condsider constructing a master ForceManager to calculate the
    % forces for all particles in the domain and then deal them between
    % objects, this will make calculation of the distances done only once
    % every loop for all objects 
    properties
        %  Names of forces acting on the particle system      
        springForce            = false;
        lennardJonesForce      = false;
        diffusionForce         = false;
        bendingElasticityForce = false;
        
        % Parameters for the forces
        bendingConst           = 0;
        springConst            = 0;
        diffusionConst         = 0;
        LJPotentialWidth       = 0;
        LJPotentialDepth       = 0;
        minParticleDistance    = 0;
        dt                     = 0; % the time to activate the force
        fixedParticleNum       = [];
                
        edges % matrices representing the edges between connected particles
        particleDistance % pairwise distance between particles
        
    end
    
    methods
        
        function obj = ForceManager(params)
            % set input parameters 
            f = fieldnames(params);
            for fIdx = 1:numel(f)
                obj.(f{fIdx}) = params.(f{fIdx});                
            end
        end
        
        function newParticlePosition = Apply(obj,particlePosition,connectivityMap,...
                                             springConst,diffusionConst,bendingConst,...
                                             LJPotentialWidth,LJPotentialDepth,...
                                             minParticleDistance,fixedParticleNum,dt)
                
            % Apply chosen forces on the N vertices/particles of an object
            % object must contain the position of its vertices in any
            % dimension
            % particlePosition should be an NxDim array of coordinates
            % particledist is the pairwise distance between the particles
            % NxN matrix
            % connectivityMap is a binary graph representing which particle
            % is connected to which other particle                        
            
            % get pair-wise distance between particles
            particleDistances    = obj.GetParticleDistance(particlePosition);
            obj.particleDistance = particleDistances;% save property
            
            % obtain the edges length
            edgeLength    = obj.GetEdgesLength(particleDistances,connectivityMap);
            
            % Apply forces
            springForces  = obj.GetSpringForce(edgeLength,springConst,connectivityMap,minParticleDistance,fixedParticleNum);
            
            % Lenard-jones force
            ljForces      = obj.GetLenardJonesForce(particlePosition,particleDistances,LJPotentialWidth,LJPotentialDepth);
            
            % Thermal (diffusion) force
            diffusionForces = obj.GetDiffusionForce(particlePosition,diffusionConst,dt);
            
            % Bending forces
            bendingForces = obj.GetBendingElasticityForce(particlePosition,connectivityMap,bendingConst);
            % The effect of applying forces is addative
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
        
        function force  = GetBendingElasticityForce(obj,particlePosition,connectivityMat,bendingConst)% needs clean up
            
            force = zeros(size(particlePosition,1),size(particlePosition,2));
            if obj.bendingElasticityForce
                % Calculate the edges vectors in all dimension
                edgeMat = GetEdgesVectors_mex(particlePosition,connectivityMat);
                force   = BendingElasticityForce_mex(edgeMat(:,:,1),...
                                                     edgeMat(:,:,2),...
                                                     edgeMat(:,:,3),...
                                                     connectivityMat,bendingConst);
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
