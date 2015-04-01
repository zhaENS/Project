classdef ForceManager<handle
    % This class collects different forces acting upon a body/particle
    % system.
    % Forces are considered to be addative, hence the order of their
    % operations is not important 
    % forceManager is a child of a physical object like a chain or a
    % particle system 
    % connectivityMap should be logical 
    % TODO: add a listener to object parameter change
    
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
                diffusionConst,springConst,LJPotentialWidth,LJPotentialDepth,minParticleDist,fixedParticleNum)
            % Apply chosen forces on the N vertices/particles of an object
            % object must contain the position of its vertices in any
            % dimension 
            % particlePosition should be an NxDim array of coordinates
            % particledist is the pairwise distance between the particles
            % NxN matrix
            % connectivityMap is a binary graph representing which particle
            % is connected to which other particle 
            
            % Calculate the edges vectors in all dimension 
             obj.edges   = obj.GetEdgesVectors(particlePosition,connectivityMap);
             
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
             bendingForces = obj.GetBendingElasticityForce(obj.edges,connectivityMap);
             % the effect of appliying forces is considered addative 
             newParticlePosition = -springForces*dt*particlePosition+...
                                      ljForces*dt+...
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
        
        function force  = GetBendingElasticityForce(obj,edgeMat,connectivityMat)% needs clean up 
            %===========    
            % Calculate the force on each particle given its connectivity 
            numParticles = size(connectivityMat,1);
            dimension    = size(edgeMat,3);
            force = zeros(numParticles,dimension);
            for pIdx = 1:numParticles
               edgeVec =  edgeMat(pIdx,connectivityMat(pIdx,:),:);
%                calculate the pairwise angles between the connected
%                particles
               n     = size(edgeVec,2);% number of edge pairs
               if n>=2
                next  = 1;
                
                numPairs = factorial(n)/(factorial(2)*factorial(n-2));
                pairs    = zeros(numPairs,2);
                forceVec  = zeros(numPairs,dimension);
                angle    = zeros(1,numPairs);
                 for a1Idx = 1:n                
                     for a2Idx = a1Idx+1:n
                         vec1          = edgeVec(1,a1Idx,:);
                         vec1          = vec1(:);
                         vec2          = edgeVec(1,a2Idx,:);
                         vec2          = vec2(:);
                         angle(next)   = acos(sum(vec1.*vec2)./(norm(vec1(:))*norm(vec2)));
                         % calculate the resulting force vector 
                         forceVec(next,:) = angle(next)./pi *( (vec1-vec2)/2);
                         pairs(next,1) = a1Idx;
                         pairs(next,2) = a2Idx;
                         next          = next+1;
                     end                     
                 end
%                  calculate the resulting force vector for the specific
%                  particle 
                force(pIdx,:) = sum(forceVec,1);
                 % calculate the force based on the pair of vectors 
               end 
            end
            
            
            
            
            
%             %%%%%%%%%%%%%
%                 % get the coordinates of the first vector comprising the
%                 % angle
%                 t1x     = obj.positions.springs.x(obj.connectionMap.indices.in.linear.bead1ToBead2);
%                 t1y     = obj.positions.springs.y(obj.connectionMap.indices.in.linear.bead1ToBead2);
%                 t1z     = obj.positions.springs.z(obj.connectionMap.indices.in.linear.bead1ToBead2);
%                 % get the coordinates of the second vector  comprising the
%                 % angle
%                 t2x     = obj.positions.springs.x(obj.connectionMap.indices.in.linear.bead3ToBead2);
%                 t2y     = obj.positions.springs.y(obj.connectionMap.indices.in.linear.bead3ToBead2);
%                 t2z     = obj.positions.springs.z(obj.connectionMap.indices.in.linear.bead3ToBead2);
%                
%                 % calculate the the force (normalized) given by the angle
%                 FValue = (pi- obj.positions.springs.angleBetweenSprings(obj.connectionMap.indices.in.linear.beadTriplets))/pi;
%                 % direct the force of each bead toward the base of the
%                 % triangle formed by the two edges (vectors) comprising the
%                 % angle
%                 tempFx  = (t1x-(0.5*(t2x-t1x))).*FValue;
%                 tempFy  = (t1y-(0.5*(t2y-t1y))).*FValue;
%                 tempFz  = (t1z-(0.5*(t2z-t1z))).*FValue;  
% %                 % ====== debug================
% %                  hold on 
% %                  line('XData',obj.positions.beads.cur.x(2),...
% %                       'YData',obj.positions.beads.cur.y(2),...
% %                       'ZData',obj.positions.beads.cur.z(2),...
% %                       'Marker','o',...
% %                       'MarkerFaceColor','r',...
% %                       'Tag','point')
% %                  % plot the vectors 
% %                  quiver3(obj.positions.beads.cur.x(2),...
% %                      obj.positions.beads.cur.y(2),...
% %                      obj.positions.beads.cur.z(2),t1x,t1y,t1z)
% %                  quiver3(obj.positions.beads.cur.x(2),...
% %                      obj.positions.beads.cur.y(2),...
% %                      obj.positions.beads.cur.z(2),t2x,t2y,t2z)
% %                   quiver3(obj.positions.beads.cur.x(2),...
% %                      obj.positions.beads.cur.y(2),...
% %                      obj.positions.beads.cur.z(2),tempFx,tempFy,tempFz)
% %                  hold off
% %                   
% %                  %=== end debug ====================
%                  
%                 % for each bead, calculate the force as the sum of forces
%                 % exerted by the angles
% %                 n       = [1:obj.params.numBeads,1:abs(numel(obj.connectionMap.beadTriplets(:,2))-obj.params.numBeads)];
% %                 a       = sparse(obj.params.numBeads,1:numel(n),0);     
%                 a       = sparse(size(obj.positions.springs.angleBetweenSprings,1),size(obj.positions.springs.angleBetweenSprings,2));
%                 ind     = obj.connectionMap.indices.uniqueBeadTripletsLinear;
%                 a       = full(a);
%                 a(ind)  = tempFx;
% 
%                 Force.x = sum(a,2);
%                 Force.x(obj.params.fixedBeadNum) = 0;
%                 a(ind)  = tempFy;
%                 Force.y = sum(a,2);
%                 Force.y(obj.params.fixedBeadNum) = 0;
%                 a(ind)  = tempFz;
%                 Force.z = sum(a,2);
%                 Force.z(obj.params.fixedBeadNum) = 0;
%                 obj.forces.bending = Force;       
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
        
        function angles = CalculateAngleBetweenSprings(springPos)% needs cleanup
                        
            linInd1 = obj.connectionMap.indices.in.linear.bead1ToBead2;
            linInd2 = obj.connectionMap.indices.in.linear.bead3ToBead2;
            c1      = full([obj.positions.springs.x(linInd1),obj.positions.springs.y(linInd1),obj.positions.springs.z(linInd1)]);
            c2      = full([obj.positions.springs.x(linInd2),obj.positions.springs.y(linInd2),obj.positions.springs.z(linInd2)]);
%             c1 = c1(sum(c1,2)~=0,:);
%             c2 = c2(sum(c2,2)~=0,:);
            sAngle  = acos(sum(c1.*c2,2)./(sqrt(sum(c1.^2,2)).*sqrt(sum(c2.^2,2))));
            linInd3 = obj.connectionMap.indices.in.linear.beadTriplets;%(obj.connectionMap.beadTriplets(:,3)-1)*obj.params.numBeads^2+linInd1;
            obj.positions.springs.angleBetweenSprings(linInd3) = sAngle;
        end
    end
    
end
