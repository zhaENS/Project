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
        
        function newParticlePosition = Apply(obj,particleDistances,particlePosition,connectivityMap,...
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
%             particleDistances    = obj.GetParticleDistance(particlePosition);
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
                    
    end
    
    methods (Static)
        
         function newParticlePosition = ApplyInternalForces(particlePosition,particleDistances,connectivityMap,...
                                                               springForce,bendingElasticityForce,...
                                                               springConst,bendingConst,...
                                                               minParticleDistance,fixedParticleNum,dt)
                                            
             % Apply object's internal forces to get the new position of
             % its vertices, represented by newParticlePosition                         
            
            % Spring force
            springForces  = ForceManager.GetSpringForce(springForce,particlePosition,particleDistances,springConst,...
                                               connectivityMap,minParticleDistance,fixedParticleNum);
            
            % Bending forces
            bendingForces = ForceManager.GetBendingElasticityForce(bendingElasticityForce,particlePosition,...
                                                                   connectivityMap,bendingConst,fixedParticleNum);
            %zero-out fixed particle position 
%             bendingForces(fixedParticleNum,:) = 0;
            dx = (springForces+ bendingForces)*dt;
            
            % zero out forces for fixed beads
%             dx(fixedParticleNum,:) = 0;
            newParticlePosition = particlePosition+dx;
                                   
         end
        
        function newParticlePosition = ApplyExternalForces(particlePosition,...
                                                            particleDistances,diffusionConst,...
                                                            ljForce,diffusionForce,...                                                            
                                                            LJPotentialWidth,LJPotentialDepth,...
                                                            fixedParticleNum,dt)
             % Apply external forces on an object to get the new position
             % for its vertices represented by newParticlePosition
             
             % Lenard-jones force             
            ljForces      = ForceManager.GetLenardJonesForce(ljForce,particlePosition,particleDistances,LJPotentialWidth,LJPotentialDepth);
            
            % Thermal (diffusion) force
            diffusionForces = ForceManager.GetDiffusionForce(diffusionForce,particlePosition,diffusionConst,dt);
            
            % zero-out forces for fixed particles
             dx = (ljForces*dt+ diffusionForces);
             dx(fixedParticleNum,:) = 0;
             
             % get new position
             newParticlePosition = particlePosition+dx;
                                                                                 
        end
        
       function newParticlePosition = ApplyCompositeInternalForces(pos,particleDistance,connectivityMap,...
                                                                        springForce,bendingElasticityForce,...
                                                                        springConst,bendingConst,...                                                         
                                                                         minParticleDistance,fixedParticleNum,dt)
                                                     
            % Apply forces on a composite structure composed of several chains with different parameters                          
            % pariclePosition is the position for each chain in cell array            
            % particleDistance is the pair-wise distance between particles
            % in object by the order of thier appearance
            % connectivityMap- is a binary matrix of connectivity between
            % all particles in the composite object
            % springForce is an N by 1 binary flag vector for each object 
            % bendingElasticityForce is an N by 1 binary flag vector for
            % each object 
            % springConst is an N by N matrix of spring constant for all
            % particles in the composite object 
            % bendingConst- is an N by 1 vector of bending constants for
            % each member of the composite object 
            % minParticleDistance is an N by N matrix for minimal distance
            % between particles in the composite object 
            % fixedParticleNum is a vector with number of particles
            % remaining fixed             
            cNb = 0;
            indsObj  = cell(1,numel(pos));
            particlePosition = [];
            for oIdx = 1:numel(pos)
             indsObj{oIdx} = (cNb+1):(cNb)+size(pos{oIdx},1);         
             particlePosition = [particlePosition;pos{oIdx}];
%              cNb = cNb+size(pos{oIdx},1);             
            end
            % Spring force
            if any(springForce)
                springForceFlag =true;
            else
                springForceFlag = false;
            end 
            
            springForces  = ForceManager.GetSpringForce(springForceFlag,particlePosition,particleDistance,...
                                            springConst,connectivityMap,minParticleDistance,fixedParticleNum);            
            % zero-out forces for object with no spring force
            % active
            springForces([indsObj{~springForce}],:)=0;
                    
            % Bending forces
             if any(bendingElasticityForce)
                bendingElasticityForceFlag = true;
            else
                bendingElasticityForceFlag = false;
             end            
            bendingForces = ForceManager.GetBendingElasticityForce(bendingElasticityForceFlag,particlePosition,connectivityMap,bendingConst(1),fixedParticleNum);
            % zero out forces for object with no bending elasticity force
            % active
            bendingForces([indsObj{~bendingElasticityForce}],:) = 0;
            
            % The effect of applying forces is addative
            dx =(springForces+bendingForces)*dt;                               
            % zero-out fixed particles position change
%             dx(fixedParticleNum,:) = 0;         
            
            newParticlePosition = particlePosition+dx; 
                                   
                                       
        end    
        
        function force  = GetSpringForce(springForce,particlePosition,particleDist,springConst,connectivityMap,minParticleDist,fixedParticleNum)
            % Calculate the spring force between N particles in any dimension M.
            % particleDist    - NxN matrix of pairwise particle distances
            % springConst     - NxN double matrix of spring constants
            % connectivityMap - NxN binary matrix which defines the connectivity between particles
            % minParticleDist - minimal distance between particles
            % fixedParticleNum - particles in the system which do not move
            force = zeros(size(particlePosition));
            if springForce
                 % obtain the edges length
%                 edgeLength    = ForceManager.GetEdgesLength(particleDist,connectivityMap);
%                 connectivityMap           = (connectivityMap);
                L                         = (1-minParticleDist./particleDist).*connectivityMap;
                L(~connectivityMap)       = 0;
                sumForces                 = sum(L,2);
                force                     = -springConst.*(diag(sumForces)-L); % set the maindiagonal                
                force                     = force*particlePosition;
                force(fixedParticleNum,:) = 0;% zero out forces for fixed particles
%                 force(:,fixedParticleNum) = 0;% zero out forces for fixed particles
            end
        end
        
        function force  = GetBendingElasticityForce(bendingElasticityForce,particlePosition,connectivityMat,bendingConst,fixedParticleNum)
            % Get bending elasticity force
            force = zeros(size(particlePosition,1),size(particlePosition,2));
            if bendingElasticityForce
                % Calculate the edges vectors in all dimension
%                  edgeMat = GetEdgesVectors(particlePosition,connectivityMat);
%                  force   = BendingElasticityForce(edgeMat(:,:,1),...
%                                                      edgeMat(:,:,2),...
%                                                      edgeMat(:,:,3),...
%                                                      connectivityMat,bendingConst);
                edgeMat = GetEdgesVectors_mex(particlePosition,connectivityMat);
                force   = BendingElasticityForce_mex(edgeMat(:,:,1),...
                                                     edgeMat(:,:,2),...
                                                     edgeMat(:,:,3),...
                                                     connectivityMat,bendingConst);
              % zero out forces for
              % fixed particles
               force(fixedParticleNum,:) = 0;
            end            
        end
        
        function force  = GetDiffusionForce(diffusionForce,particlePosition,diffusionConst,dt)
            % get thermal (diffusion) force
            force = zeros(size(particlePosition));
            if diffusionForce
                force = randn(size(particlePosition))*sqrt(2*diffusionConst*dt);
            end            
        end
        
        function force  = GetLenardJonesForce(ljForce,particlePosition,particleDist,LJPotentialWidth,LJPotentialDepth)
            % calculate Lenard jones force between particles
            force = zeros(size(particlePosition));
            if ljForce
%                 force = LennardJones(particlePosition,particleDist,LJPotentialWidth,LJPotentialDepth);
                force = LennardJones_mex(particlePosition,particleDist,LJPotentialWidth,LJPotentialDepth);
            end
        end
        
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
%             particleDistance = pdist2(particlePosition,particlePosition);
            particleDistance = pdist2mex(particlePosition',...
                particlePosition','euc',[],[],[]);
        end
        
    end
    
end
