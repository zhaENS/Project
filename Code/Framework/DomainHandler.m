classdef DomainHandler<handle
    % A domain class to be used in the polymer simulation 
    
    % TODO: allow insertion of a polygonal domain or a mesh from meshlab in
    % a *.stl file format, expand Reflect function accordingly 
    % TODO: fix reflection for all domain shapes
    % TODO: prepare reflection function such that it will output location on
    %      the boundary and allow attachment
    % TODO: complete reflection 'out' in all domain shapes
    % TODO: add polygon to the list of domains 
    properties 
        params
        points
        numDomains
    end    
    
    methods 
        function obj = DomainHandler(domainParams)
             % class constructor 
            obj.params     = domainParams;
            obj.numDomains = numel(domainParams);                      
        end                                
        
        function SetDefaultParams(obj)%obsolete
            obj.params.domainShape    = 'Sphere'; % sphere| cylinder | twoPlates | open
            obj.params.domainWidth    = 1;
            obj.params.domainHeight   = 50;
            obj.params.dimension      = 3;
            obj.params.reflectionType = 'preserveEnergy'; 
            obj.params.domainCenter   = [0 0 0];
            obj.params.parentAxes     = [];                        
        end
        
        function SetInputParams(obj,varargin)
            if mod(numel(varargin{:}),2)~=0
                error('input must come in name-value pairs')
            else 
                for vIdx = 1:2:numel(varargin)
                    obj.params.(varargin{vIdx}) = varargin{vIdx+1};
                end
            end
            
            obj.handles.graphical.mainAxes = obj.params.parentAxes;
        end        
                
        function newParticlePosition = Step(obj, prevParticlePosition,curParticlePosition,particleDist,fixedParticleNum,particlesOnBoundary,domainInds,dt)                                                
               % Apply forces on particles in the domain to obtain their
               % new position 
               % currently this function supports domains which do not
               % overlap or contained domains in which one is completely
               % permeable 
               
               % experimental- apply forces iteratively from each domain 
               for dIdx = 1:obj.numDomains
                 
                   domainNumber = dIdx;
                   p    = domainInds(particlesOnBoundary);
                   pInd = (p==dIdx);%obtain the index of particlesOnBoundary for each domain;
                   dp           = obj.params(dIdx);
                   fp           = dp.forceParams;
                   dInds = (domainInds==dIdx);
                   %if has more than one domain,as the particlesOnBoundary
                   %is a global number,change the index to fit current
                   %domain;
%                    if dIdx >1
%                    particlesOnBoundary(pInd) = particlesOnBoundary(pInd)-numel(find(dInds~=0));
%                    end                  
                   cp = ForceManager.ApplyExternalForces(curParticlePosition,particleDist,...
                                                    fp.diffusionConst,fp.lennardJonesForce,...
                                                    fp.diffusionForce,fp.morseForce,fp.mechanicalForce,...
                                                    fp.LJPotentialWidth,fp.LJPotentialDepth,...
                                                    fp.morsePotentialDepth,fp.morsePotentialWidth,fp.morseForceType,...
                                                    fp.mechanicalForceCenter,fp.mechanicalForceDirection,fp.mechanicalForceMagnitude,...
                                                    fp.minParticleEqDistance,fixedParticleNum,...
                                                    particlesOnBoundary(pInd),dt);
                    curParticlePosition(dInds,:)= cp(dInds,:);
                                                      
                    
                   if ~strcmpi(dp.reflectionType,'off')% if reflection is set to on 
                     [~,curParticlePosition(dInds,:)] = obj.Reflect(prevParticlePosition(dInds ,:),curParticlePosition(dInds,:),domainNumber);
                   end
               end
               
               newParticlePosition = curParticlePosition;
               
        end
                     
        function [prevPos,curPos,inFlag] = Reflect(obj,pos1,pos2,domainNumber,reflectDirection)%TODO: complete reflection out in all domain shapes
              % Reflect a particle previously at pos1 and currently at
              % pos2 depending on the domain shape 
              % pos1 and pos2  are  NxDim arrays of particle positions 
              % domainNumber is the serial number of the domain registered in the class
              % reflectDirection is either 'in' or 'out' indicating the direction of reflection,
              % inside the domain or outside the domain
              
              if ~exist('reflectDirection','var')
                  reflectDirection = 'in'; % assume reflection in the domain by default
              end
              
                curPos  = pos2;
                prevPos = pos1;
                
                numBeads = size(pos1,1);
                inFlag   = true(numBeads,1);         
                                                 
                [inIdx,onIdx] = obj.InDomain(pos2,domainNumber); % find all particles in the domain 
                outIdx        = find(~(inIdx|onIdx));       % find all particles outside the domain 
%                 dc       = obj.params(domainNumber).domainCenter;
                % For each particle outside the domain, reflect until it is
                % inside the domain 
            for oIdx = 1:numel(outIdx)
                 count          = 1; % counter for reflections 
                 ind            = outIdx(oIdx);
                 curPos(ind,:)  = pos2(ind,:);
                 prevPos(ind,:) = pos1(ind,:);
                while ~obj.InDomain(curPos(ind,:),domainNumber)
                    % Reflect the particle
                    if strcmpi(obj.params(domainNumber).domainShape,'open')
                        % Do nothing
                    elseif strcmpi(obj.params(domainNumber).domainShape,'sphere')
                        
                        % Find intersection with the domain 
                        intersectionPoint = obj.FindIntersectionPoint(prevPos(ind,:),curPos(ind,:),domainNumber);
                        
                        if ~obj.InDomain(prevPos(ind,:),domainNumber)
                            error('something is wrong with the intersection point function')
                        end
                        
                        if ~isempty(intersectionPoint)
                            domainNorm  = obj.GetDomainNormal(intersectionPoint,domainNumber,reflectDirection);
                            
                            % Calculate reflected ray
                            pp  = prevPos(ind,:);
                            cp  = curPos(ind,:);
                            di  = (intersectionPoint-pp)/sqrt(sum(intersectionPoint-pp).^2);
                            ds  = -(2*dot(domainNorm,di)*domainNorm-di);
                            t   = sqrt(sum(cp-intersectionPoint).^2);
                            n   = intersectionPoint+ t*ds;% the new position
                            
                            biasNorm = (n-intersectionPoint)/norm(n-intersectionPoint);

                                curPos(ind,:)  = n;
                                % To avoid numerical error, move the prev point
                                % slightly on the vector between the new point and the
                                % intersectionPoint
                                prevPos(ind,:) = intersectionPoint+(obj.params(domainNumber).domainWidth/1e5)*biasNorm;
                                
%                             end
                        else
                            disp('no intersection point')
                        end
                        
                    elseif strcmpi(obj.params.domainShape,'cylinder');
                        % find the intersection in 2D first
                        intersectionPoint = obj.FindIntersectionPoint([prevPos(:,1) prevPos(:,2) zeros(size(prevPos,1),1)],[curPos(:,1) curPos(:,2) zeros(size(curPos,1),1)],domainNumber);
                        
                        if ~isempty(intersectionPoint)
                            theta  = subspace([curPos(:,1)-prevPos(:,1) curPos(:,2)-prevPos(:,2) curPos(:,3)-prevPos(:,3)]',[intersectionPoint(1)-prevPos(:,1), intersectionPoint(2)-prevPos(:,2),prevPos.z]');
                            intersectionPoint(end) = obj.params(domainNumber).domainWidth/cos(theta);% z(prevPos.x,curPos.x,intersectionPoint(1),prevPos.y,curPos.y,intersectionPoint(2),prevPos.z,curPos.z);
                            domainNorm  = obj.GetDomainNormal(intersectionPoint,domainNumber);
                            cp      = [curPos.x curPos.y curPos.z];
                            pp      = [prevPos.x prevPos.y prevPos.z];
                            d       = (cp-pp)/sqrt(sum(cp-pp).^2);
                            r       = d-2*dot(d,domainNorm)*domainNorm;
                            n       = intersectionPoint+r*sqrt(sum(intersectionPoint-cp).^2)/sqrt(sum(r.^2));
                            curPos  = n;
                            prevPos = intersectionPoint;                            
                        end
                        
                    elseif strcmpi(obj.params(domainNumber).domainShape,'twoPlates')
                        % The two plates are arranged vertically in parrallel
                        % the distance between the plates is
                        % obj.params.domainWidth
                        intersectionPoint = obj.FindIntersectionPoint(prevPos,curPos,domainNumber);% the intersection point is made sure to lay on the boundary
                        if ~isempty(intersectionPoint)
                            domainNorm        = obj.GetDomainNormal(intersectionPoint,domainNumber);
                            d                 = (curPos-prevPos)/sqrt(sum(curPos-prevPos).^2);
                            r                 = d-2*dot(d,domainNorm)*domainNorm;
                            n                 = intersectionPoint+r*sqrt(sum(intersectionPoint-curPos).^2)/sqrt(sum(r.^2));
%                             obj.PlotReflection(prevPos,curPos,intersectionPoint,n,domainNorm);
                            curPos  = n;
                            prevPos = intersectionPoint;
                        end
                        
                    else
                        prevPos = pp;
                    end
                    
                    count = count+1;
                    if count>obj.params(domainNumber).maxReflectionsPerParticle
                        error('%s%s%s','Reflection count is bigger than ',num2str(obj.params(domainNumber).maxReflectionsPerParticle),...
                            '- Terminating');                        
                    end
                end
            end
        end
        
        function intersectionPoint = FindIntersectionPoint(obj,prevPos,curPos,domainNumber)

            if strcmpi(obj.params(domainNumber).domainShape,'sphere') || strcmpi(obj.params(domainNumber).domainShape,'cylinder')
                    % Find the intersection of a point and the surface of the domain
                     dc = obj.params(domainNumber).domainCenter;
                     R  = obj.params(domainNumber).domainWidth;% domain radius                    
                     A  = prevPos-dc;
                     B  = curPos-dc; 
                     C  = (B-A);%./norm(B-A);
                     
                     gamma = dot(A,A);
                     alpha = dot (A,B)-gamma;
                     beta  = dot(C,C);
                     t(1)  = (-alpha+sqrt(alpha^2 -beta*(gamma-R^2)))/(beta);
                     t(2)  = (-alpha-sqrt(alpha^2 -beta*(gamma-R^2)))/(beta);
                     
                     
                      % === end test ====
                     % Take only the positive root smaller than 1 
                    t    = t(t>1e-13&t<=1); 
                    if ~isempty(t) && isreal(t)   
                        if numel(t)>1
                            disp('two roots')
                            
                        end
                       intersectionPoint = prevPos+min(t)*(C);
%                        intersectionPoint = intersectionPoint+dc;
                    else
                        intersectionPoint = [];
                    end
            elseif strcmpi(obj.params(domainNumber).domainShape,'twoPlates')
                    r  = obj.params(domainNumber).domainWidth;                                    
                    c  = cross([0 1 0],[0 0 1]);
                    t1 = dot([r 0 0]-prevPos,c)/(dot(curPos-prevPos,c));
                    t2 = dot([-r 0 0] -prevPos,c)/(dot(curPos-prevPos,c));                    
                    if t1>0
                        intersectionPoint = prevPos+t1*(curPos-prevPos);
                    elseif t2>0
                        intersectionPoint = prevPos+t2*(curPos-prevPos);
                    else 
                        intersectionPoint = [];
                    end
            end
        end
        
        function [inIdx, onIdx] = InDomain(obj,vecIn,domainNumber)
            % Check if points in vecIn are inside domain(domainNumber)            
            % the output is a binary vector with 1 if the point is inside
            % the domain, 0 otherwise.           
            if ~exist('domainNumber','var')
                domainNumber= 1;
            end
            
            numParticles = size(vecIn,1);
%             for dIdx = 1:obj.numDomains
                inIdx  = true(numParticles,1);
                onIdx  = false(numParticles,1);
                if strcmpi(obj.params(domainNumber).domainShape,'Sphere')
                    % the vector norm
                    dc    = obj.params(domainNumber).domainCenter;
                    n     = sqrt(sum(bsxfun(@minus,vecIn,dc).^2,2));
                    onIdx = ((n-obj.params(domainNumber).domainWidth).^2)<eps;
                    inIdx = (n+1e-14<(obj.params(domainNumber).domainWidth));
                   if any(onIdx& inIdx)
                       disp('stop')
                   end
                elseif strcmpi(obj.params(domainNumber).domainShape,'cylinder')                
                    % the vector norm
                    dc    = obj.params(domainNumber).domainCenter;
                    dc    = [dc(1) dc(2) 0];
                    n     = sqrt(sum(bsxfun(@minus,vecIn(:,1:2),dc(:,1:2)).^2,2));
                    inIdx = n<=obj.params(domainNumber).domainWidth;
                elseif strcmpi(obj.params(domainNumber).domainShape,'twoPlates')
                    inIdx = vecIn.x<obj.params(domainNumber).domainWidth;
                elseif strcmpi(obj.params(domainNumber).domainShape,'open')
                    % Do nothing
                end            
%             end
        end
        
        function domainNorm = GetDomainNormal(obj,point,domainNumber,normDirection)%TODO: finish inserting domain direction to TwoPlates
           % get the normal vector to of the domain at location indicated by point
           % the domain number is the serial number of the domain registered in teh class. 
           % normDirection is either 'in', or 'out', indicating the normals pointing into the domain
           % or outside            
           % The normal vectors' norms are normalized to unity
           
           if ~exist('domainNumber','var')
               domainNumber = 1;
           end
            if exist('normDirection','var')
                normDirection = 'in';
            end
            
            if strcmpi(normDirection,'in')
                direction = 1;
            else
                direction = -1;
            end
            
           if strcmpi(obj.params(domainNumber).domainShape,'sphere')
               domainNorm = direction* obj.params(domainNumber).domainCenter-point;
               domainNorm = domainNorm/sqrt(sum(domainNorm.^2));
               
           elseif strcmpi(obj.params(domainNumber).domainShape,'cylinder')
               domainNorm    = direction*(obj.params(domainNumber).domainCenter-point);
               domainNorm(3) = 0;
               domainNorm = domainNorm/sqrt(sum(domainNorm.^2));
           elseif strcmpi(obj.params(domainNumber).domainShape,'twoPlates')
               if point(1)< obj.params(domainNumber).domainWidth+1e-7 && point(1)> obj.params(domainNumber).domainWidth-1e-7 
                   domainNorm = -[1,0 0];
               elseif point(1)<-obj.params(domainNumber).domainWidth+1e-7 && point(1)>-obj.params(domainNumber).domainWidth-1e-7
                   domainNorm = [1 0 0];
               end
           elseif  strcmpi(obj.params(domainNumber).domainShape,'open')
               domainNorm = [];
           end
          
           
        end
        
        function curParticlesPositions = MoveDomain(obj,curPos,domainNumber)
             if ~exist('domainNumber','var')
               domainNumber = 1;
             end
             
             if strcmpi(obj.params(domainNumber).moveDomainType,'none')
                curParticlesPositions = curPos;
             else
                 if strcmpi(obj.params(domainNumber).moveDomainType,'rotate')
                 
                 end
                 if strcmpi(obj.params(domainNumber).moveDomainType,'straight')
                 
                 
                 end
                 if strcmpi(obj.params(domainNumber).moveDomainType,'rotate&straight')
                 
                 
                 end
             end
        end
        
        
        
        function boundaryPoints = GetRandomBoundarySample(obj,numPoints,domainNumber)%TODO: add support for other domain shapes
            % numPoints random sample of the domain 
            if ~exist('domainNumber','var')
                domainNumber = 1;
            end
            
            boundaryPoints = [];
            if strcmpi(obj.params(domainNumber).domainShape,'open')
                % normaly random distributed points
                boundaryPoints = randn(numPoints,obj.params(domainNumber).dimension);
            elseif strcmpi(obj.params(domainNumber).domainShape,'sphere')
                % There is an analytical representation of the sphere 
                % random sample 
                u     = rand(numPoints,1);
                v     = rand(numPoints,1);
                phi   = 2*pi*u;
                theta = acos(2*v -1);
                boundaryPoints(:,1) = obj.params(domainNumber).domainWidth.*sin(theta).*cos(phi)+obj.params(domainNumber).domainCenter(1);
                boundaryPoints(:,2) = obj.params(domainNumber).domainWidth.*sin(theta).*sin(phi)+obj.params(domainNumber).domainCenter(2);
                boundaryPoints(:,3) = obj.params(domainNumber).domainWidth.*cos(theta)+obj.params(domainNumber).domainCenter(3);
                
            end
            
        end
        
    end
    
end