classdef DomainHandler<handle
    % A domain class to be used in the polymer simulation 
    
    % TODO: allow insertion of a polygonal domain or a mesh from meshlab in
    % a *.stl file format, expand Reflect function accordingly 
    % TODO: fix reflection for all domain shapes
    % TODO: prepare reflection function such that it will output location on
    %      the boundary and allow attachment
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
                
        function newParticlePosition = Step(obj, prevParticlePosition,curParticlePosition,particleDist,fixedParticleNum,dt)                                                
               % Apply forces on particles in the domain to obtain their
               % new position 
               % currently this function supports domains which do not
               % overlap or contained domains in which one is completely
               % permeable 
               
               % experimental- apply forces iteratively from each domain 
               for dIdx = 1:obj.numDomains
                   domainNumber = dIdx;
                   dp = obj.params(dIdx);
                   fp = dp.forceParams;
                   curParticlePosition = ForceManager.ApplyExternalForces(curParticlePosition,particleDist,...
                                    fp.diffusionConst,fp.lennardJonesForce,fp.diffusionForce, fp.morseForce,...
                                    fp.LJPotentialWidth,fp.LJPotentialDepth,...
                                    fp.morsePotentialDepth, fp.morsePotentialWidth,fp.morseForceType,...
                                    fp.minParticleEqDistance,fixedParticleNum,dt);
                                                                                               
                   if ~strcmpi(dp.reflectionType,'off')% if reflection is set to on 
                     [~,curParticlePosition] = obj.Reflect(prevParticlePosition,curParticlePosition,domainNumber);
                   end
               end
               
               newParticlePosition = curParticlePosition;
               
        end
             
        
        function [prevPos,curPos,inFlag] = Reflect(obj,pos1,pos2,domainNumber)
              % Reflect a particle previously at pos1 and currently at
              % pos2 depending on the domain shape 
              % pos1 and pos2  are  NxDim arrays of particle positions 
                curPos  = pos2;
                prevPos = pos1;
                
                numBeads = size(pos1,1);
                inFlag   = true(numBeads,1);         
                                     
            
                inIdx    = obj.InDomain(pos2,domainNumber); % find all particles in the domain 
                outIdx   = find(~inIdx);       % find all particles outside the domain 
                
                % For each particle outside the domain, reflect until it is
                % inside the domain 
            for oIdx = 1:numel(outIdx)
                 count    = 1; % counter for reflections 
                 ind      = outIdx(oIdx);
                 curPos(ind,:)  = pos2(ind,:);
                 prevPos(ind,:) = pos1(ind,:);
                while ~obj.InDomain(curPos(ind,:),domainNumber)
                    % Reflect the particle
                    if strcmpi(obj.params(domainNumber).domainShape,'open')
                        % Do nothing
                    elseif strcmpi(obj.params(domainNumber).domainShape,'sphere')
                        
                        % Find intersection with the domain 
                        intersectionPoint = obj.FindIntersectionPoint(prevPos(ind,:),curPos(ind,:),domainNumber);
                        
                        if all([obj.InDomain(prevPos(ind,:),domainNumber),~obj.InDomain(curPos(ind,:),domainNumber),isempty(intersectionPoint)]) %%|| sqrt(sum(intersectionPoint.^2))~=obj.params.domainWidth)
                            error('something is wrong with the intersection point function')
                        end
                        
                        if ~isempty(intersectionPoint)
                            domainNorm  = obj.GetDomainNormal(intersectionPoint,domainNumber);
                            
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
                    
                     R  = obj.params(domainNumber).domainWidth;% domain radius                    
                     A  = prevPos;
                     B  = curPos; 
                     C  = (B-A);%./norm(B-A);
                     
                     gamma = dot(A,A);
                     alpha = dot (A,B) -gamma;
                     beta  = dot(C,C);
                     t(1)= (-alpha+sqrt(alpha^2 -beta*(gamma-R^2)))/(beta);
                     t(2)= (-alpha-sqrt(alpha^2 -beta*(gamma-R^2)))/(beta);
                     
                     
                      % === end test ====
                     % Take only the positive root smaller than 1 
                    t    = t(t>0&t<=1); 
                    if ~isempty(t) && isreal(t)   
                        if numel(t)>1
                            disp('two roots')
                            
                        end
                       intersectionPoint = prevPos+min(t)*(C);
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
        
        function inIdx = InDomain(obj,vecIn,domainNumber)
            % Check if points in vecIn are inside domain(domainNumber)            
            % the output is a binary vector with 1 if the point is inside
            % the domain, 0 otherwise.           
            if ~exist('domainNumber','var')
                domainNumber= 1;
            end
            
            numParticles = size(vecIn,1);
            for dIdx = 1:obj.numDomains
                inIdx        = true(numParticles,1);
                if strcmpi(obj.params(domainNumber).domainShape,'Sphere')
                    % the vector norm
                    n     = sqrt(sum(vecIn.^2,2));
                    inIdx = (n<=(obj.params(domainNumber).domainWidth+eps));

                elseif strcmpi(obj.params(domainNumber).domainShape,'cylinder')                
                    % the vector norm
                    n     = sqrt(sum(vecIn.^2,2));
                    inIdx = n<=obj.params(domainNumber).domainWidth;
                elseif strcmpi(obj.params(domainNumber).domainShape,'twoPlates')
                    inIdx = vecIn.x<obj.params(domainNumber).domainWidth;
                elseif strcmpi(obj.params(domainNumber).domainShape,'open')
                    % Do nothing
                end            
            end
        end
        
        function domainNorm = GetDomainNormal(obj,intPoint,domainNumber)
           % Assuming the domain is a shpere 
           if ~exist('domainNumber','var')
               domainNumber = 1;
           end
           
            % the normals are facing inward
           if strcmpi(obj.params(domainNumber).domainShape,'sphere')
               domainNorm = obj.params(domainNumber).domainCenter-intPoint;
               domainNorm = domainNorm/sqrt(sum(domainNorm.^2));
               
           elseif strcmpi(obj.params(domainNumber).domainShape,'cylinder')
               domainNorm    = obj.params(domainNumber).domainCenter-intPoint;
               domainNorm(3) = 0;
               domainNorm = domainNorm/sqrt(sum(domainNorm.^2));
           elseif strcmpi(obj.params(domainNumber).domainShape,'twoPlates')
               if intPoint(1)< obj.params(domainNumber).domainWidth+1e-7 && intPoint(1)> obj.params(domainNumber).domainWidth-1e-7 
                   domainNorm = -[1,0 0];
               elseif intPoint(1)<-obj.params(domainNumber).domainWidth+1e-7 && intPoint(1)>-obj.params(domainNumber).domainWidth-1e-7
                   domainNorm = [1 0 0];
               end
           elseif  strcmpi(obj.params(domainNumber).domainShape,'open')
               domainNorm = [];
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
                boundaryPoints(:,1) = obj.params(domainNumber).domainWidth.*sin(theta).*cos(phi);
                boundaryPoints(:,2) = obj.params(domainNumber).domainWidth.*sin(theta).*sin(phi);
                boundaryPoints(:,3) = obj.params(domainNumber).domainWidth.*cos(theta);
                
            end
            
        end
        
%         function PlotReflection(obj,prevPos,curPos,intPos,newPos,domainNorm)% obsolete
%             % the new pos is the point after reflection 
%             
%             % drw a line between the prev (in) point and cur (out) point
%             line('XData',[prevPos(1) curPos(1)],...
%                  'YData',[prevPos(2) curPos(2)],...
%                  'ZData',[prevPos(3) curPos(3)],...
%                  'Color','y',...
%                  'Parent',obj.handles.graphical.mainAxes);
%              
%              % show intersection point
%    ApplyForces         line('XData',intPos(1),...
%                  'YData',intPos(2),...
%                  'ZData',intPos(3),...
%                  'Marker','o',...
%                  'MarkerEdgeColor','b',...
%                  'lineStyle','none',...
%                  'Color','m',...
%                  'Parent',obj.handles.graphical.mainAxes);
%              
%              % circle the cur point
%              line('XData',curPos(1),...
%                  'YData',curPos(2),...
%                  'ZData',curPos(3),...
%                  'Marker','o',...
%                  'MarkerEdgeColor','r',...
%                  'lineStyle','none',...
%                  'Parent',obj.handles.graphical.mainAxes);
%              
%              % circel prev point
%              line('XData',prevPos(1),...
%                  'YData',prevPos(2),...
%                  'ZData',prevPos(3),...
%                  'Marker','o',...
%                  'MarkerEdgeColor','g',...
%                  'lineStyle','none',...
%                  'Parent',obj.handles.graphical.mainAxes);
%              
%              % draw a line between the reflection point and the new (in)
%              % point
%                line('XData',[intPos(1) newPos(1)],...
%                  'YData',[intPos(2) newPos(2)],...
%                  'ZData',[intPos(3) newPos(3)],...
%                  'Marker','o',...                 
%                  'MarkerEdgeColor','b',...
%                  'lineStyle','-',...
%                  'Color','c',...
%                  'Parent',obj.handles.graphical.mainAxes);
%              
%              % show the domain normal
%              hold on 
%              quiver3(intPos(1),intPos(2),intPos(3),domainNorm(1), domainNorm(2), domainNorm(3),'Color','r')
% %              quiver3(intPos(1),intPos(2),intPos(3),newPos(1), newPos(2), newPos(3),'Color','g')
%              hold off
%         end
    end
    
end