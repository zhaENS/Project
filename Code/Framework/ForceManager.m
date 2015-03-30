classdef ForceManager<handle
    % This class collects different forces acting upon a body/particle
    % system.
    % Forces are considered to be addative, hence the order of their
    % operations is not important 
    properties
        % The name of forces acting on a particle system
        springForce
        LenardJones
        Brownian
        bending
    end
    
    methods
        function obj = ForceManager
        end
        
        function Apply(obj,object)
            % apply chosen forces on the vertices/ particles of an object
            % object must contain the position of its vertices in any
            % dimension 
            
            
        end
    end
    
    methods (Static)
        
        function force  = GetSpringForce(particleDist,springConst,connectivityMap,minParticleDist,fixedParticleNum)
            % Calculate the spring force between N particles in any dimension M.
            % particleDist    - NxN matrix of pairwise particle distances
            % springConst     - NxN double matrix of spring constants
            % connectivityMap - NxN binary matrix which defines the connectivity between particles
            % minParticleDist - minimal distance between particles
            % fixedParticleNum - particles in the system which do not move
            
            connectivityMap          = double(connectivityMap);
            % force              = springConst.*particleDist.*connectivityMap;
            particleDist = particleDist+diag(Inf*ones(1,size(particleDist,1)));
            L                         = (1-minParticleDist./particleDist).*connectivityMap;
            sumForces                 = sum(L,2);
            force                     = springConst.*(diag(sumForces)-L); % set the maindiagonal
            force(fixedParticleNum,:) = 0;% zero out forces for fixed particles
            force(:,fixedParticleNum) = 0;% zero out forces for fixed particles
        end
        
        function force  = GetLenardJonesForce(particlePosition,particleDist,LJPotentialWidth,LJPotentialDepth)% move to mex
            % calculate Lenard jones force between particles 
            
                numParticles = size(particlePosition,1);
                dimension    = size(particlePosition,2);
                force        = zeros(numParticles,dimension);
                beadsDistMat = particleDist; 
                sig          = LJPotentialWidth;
                epsilon      = LJPotentialDepth;
                bdmInv       = (beadsDistMat).^(-1);            % one over the bead distance matrix 
                d            = MatIntPower(bdmInv, 6); % matrix integer power (form Utils)
                t            = (sig^6).*d;                
                forceValue   = 24*(epsilon*bdmInv).*(-2*t.*t +t); % derivative of LJ function 
                
                forceValue(isnan(forceValue))= 0;
                for dIdx = 1:size(particlePosition,2)
                    % replicate the position vector
                    A    = particlePosition(:,dIdx);
                    siz  = numParticles;
                    B1   = A(:,ones(siz(2),1));
                    siz  = [numParticles,1];
                    A    = A';
                    B2   = A(ones(siz(1),1),:);
                    % Subtract positions to get the direction vectors
                    fd  = B1-B2;

                    force(:,dIdx) = (sum(fd.*forceValue,1)');
%                     force(:,dIdx)(obj.params.fixedBeadNum) = 0; 
                end
        
        end
        
        function force  = GetBendingElasticity()% needs clean up 
                
                % get the coordinates of the first vector comprising the
                % angle
                t1x     = obj.positions.springs.x(obj.connectionMap.indices.in.linear.bead1ToBead2);
                t1y     = obj.positions.springs.y(obj.connectionMap.indices.in.linear.bead1ToBead2);
                t1z     = obj.positions.springs.z(obj.connectionMap.indices.in.linear.bead1ToBead2);
                % get the coordinates of the second vector  comprising the
                % angle
                t2x     = obj.positions.springs.x(obj.connectionMap.indices.in.linear.bead3ToBead2);
                t2y     = obj.positions.springs.y(obj.connectionMap.indices.in.linear.bead3ToBead2);
                t2z     = obj.positions.springs.z(obj.connectionMap.indices.in.linear.bead3ToBead2);
               
                % calculate the the force (normalized) given by the angle
                FValue = (pi- obj.positions.springs.angleBetweenSprings(obj.connectionMap.indices.in.linear.beadTriplets))/pi;
                % direct the force of each bead toward the base of the
                % triangle formed by the two edges (vectors) comprising the
                % angle
                tempFx  = (t1x-(0.5*(t2x-t1x))).*FValue;
                tempFy  = (t1y-(0.5*(t2y-t1y))).*FValue;
                tempFz  = (t1z-(0.5*(t2z-t1z))).*FValue;  
%                 % ====== debug================
%                  hold on 
%                  line('XData',obj.positions.beads.cur.x(2),...
%                       'YData',obj.positions.beads.cur.y(2),...
%                       'ZData',obj.positions.beads.cur.z(2),...
%                       'Marker','o',...
%                       'MarkerFaceColor','r',...
%                       'Tag','point')
%                  % plot the vectors 
%                  quiver3(obj.positions.beads.cur.x(2),...
%                      obj.positions.beads.cur.y(2),...
%                      obj.positions.beads.cur.z(2),t1x,t1y,t1z)
%                  quiver3(obj.positions.beads.cur.x(2),...
%                      obj.positions.beads.cur.y(2),...
%                      obj.positions.beads.cur.z(2),t2x,t2y,t2z)
%                   quiver3(obj.positions.beads.cur.x(2),...
%                      obj.positions.beads.cur.y(2),...
%                      obj.positions.beads.cur.z(2),tempFx,tempFy,tempFz)
%                  hold off
%                   
%                  %=== end debug ====================
                 
                % for each bead, calculate the force as the sum of forces
                % exerted by the angles
%                 n       = [1:obj.params.numBeads,1:abs(numel(obj.connectionMap.beadTriplets(:,2))-obj.params.numBeads)];
%                 a       = sparse(obj.params.numBeads,1:numel(n),0);     
                a       = sparse(size(obj.positions.springs.angleBetweenSprings,1),size(obj.positions.springs.angleBetweenSprings,2));
                ind     = obj.connectionMap.indices.uniqueBeadTripletsLinear;
                a       = full(a);
                a(ind)  = tempFx;

                Force.x = sum(a,2);
                Force.x(obj.params.fixedBeadNum) = 0;
                a(ind)  = tempFy;
                Force.y = sum(a,2);
                Force.y(obj.params.fixedBeadNum) = 0;
                a(ind)  = tempFz;
                Force.z = sum(a,2);
                Force.z(obj.params.fixedBeadNum) = 0;
                obj.forces.bending = Force;       
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
