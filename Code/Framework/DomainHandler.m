classdef DomainHandler<handle
    % A domain class to be used in the polymer simulation 
    
    % TODO: allow insertion of a polygonal domain or a mesh from meshlab in
    % a *.stl file format, expand Reflect function accordingly 
    % TODO: fix reflection for all domain shapes
    properties 
        handles
        params
        points
    end    
    
    methods 
        function obj = DomainHandler(domainParams)
             % class constructor 
            obj.params = domainParams;
                      
            obj.handles.graphical.domain   = [];
            obj.handles.graphical.light    = [];  
            obj.handles.graphical.mainAxes = obj.params.parentAxes;            
        end                        
        
        
        function SetDefaultParams(obj)%obsolete
            obj.params.domainShape    = 'Sphere'; % sphere| cylinder | twoPlates | none
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
        
        function newParticlePosition = ApplylForces(obj,prevParticlePosition,...
                                                        curParticlePosition,...
                                                        particleDist,...
                                                        ljForce,diffusionForce,...
                                                        diffusionConst,LJPotentialWidth,LJPotentialDepth,...
                                                        fixedParticleNum,dt)
                                                    
            % Apply forces on object in the domain, reflect them if
            % neccessary and return their new position 
            newParticlePosition = ForceManager.ApplyExternalForces(curParticlePosition,particleDist,diffusionConst,...
                                                                  ljForce,diffusionForce,...
                                                                  LJPotentialWidth,LJPotentialDepth,...
                                                                  fixedParticleNum,dt);
                                 
            [~,newParticlePosition] = obj.Reflect(prevParticlePosition,newParticlePosition);
                           
        end        
        
        function [prevPos,curPos,inFlag] = Reflect(obj,pos1,pos2)
              % Reflect a particle previously at pos1 and currently at
              % pos2 depending on the domain shape 
              % pos1 and pos2  are  NxDim arrays of particle positions 
                curPos  = pos2;
                prevPos = pos1;
                
                numBeads = size(pos1,1);
                inFlag   = true(numBeads,1);         
                                     
            
                inIdx    = obj.InDomain(pos2); % find all particles in the domain 
                outIdx   = find(~inIdx);       % find all particles outside the domain 
                
                % For each particle outside the domain, reflect until it is
                % inside the domain 
            for oIdx = 1:numel(outIdx)
                 count    = 1; % counter for reflections 
                 curPos(outIdx(oIdx),:)  = pos2(outIdx(oIdx),:);
                 prevPos(outIdx(oIdx),:) = pos1(outIdx(oIdx),:);
                while ~obj.InDomain(curPos(outIdx(oIdx),:))
                    % Reflect the particle
                    if strcmpi(obj.params.domainShape,'none')% change 'none' to 'open'
                        % Do nothing
                    elseif strcmpi(obj.params.domainShape,'sphere')
                        
                        % Find intersection with the domain 
                        intersectionPoint = obj.FindIntersectionPoint(prevPos(outIdx(oIdx),:),curPos(outIdx(oIdx),:));
                        
                        if all([obj.InDomain(prevPos(outIdx(oIdx),:)),~obj.InDomain(curPos(outIdx(oIdx),:)),isempty(intersectionPoint)]) %%|| sqrt(sum(intersectionPoint.^2))~=obj.params.domainWidth)
                            error('something is wrong with the intersection point function')
                        end
                        
                        if ~isempty(intersectionPoint)
                            domainNorm  = obj.GetDomainNormal(intersectionPoint);
                            
                            % Calculate reflected ray
                            pp  = prevPos(outIdx(oIdx),:);
                            cp  = curPos(outIdx(oIdx),:);
                            di  = (intersectionPoint-pp)/sqrt(sum(intersectionPoint-pp).^2);
                            ds  = -(2*dot(domainNorm,di)*domainNorm-di);
                            t   = sqrt(sum(cp-intersectionPoint).^2);
                            n   = intersectionPoint+ t*ds;% the new position
                            
                            biasNorm = (n-intersectionPoint)/norm(n-intersectionPoint);

                                curPos(outIdx(oIdx),:)  = n;
                                % To avoid numerical error, move the prev point
                                % slightly on the vector between the new point and the
                                % intersectionPoint
                                prevPos(outIdx(oIdx),:) = intersectionPoint+(obj.params.domainWidth/1e10)*biasNorm;
                                
%                             end
                        else
                            disp('no intersection point')
                        end
                        
                    elseif strcmpi(obj.params.domainShape,'cylinder');
                        % find the intersection in 2D first
                        intersectionPoint = obj.FindIntersectionPoint([prevPos(:,1) prevPos(:,2) zeros(size(prevPos,1),1)],[curPos(:,1) curPos(:,2) zeros(size(curPos,1),1)]);
                        
                        if ~isempty(intersectionPoint)
                            theta =subspace([curPos(:,1)-prevPos(:,1) curPos(:,2)-prevPos(:,2) curPos(:,3)-prevPos(:,3)]',[intersectionPoint(1)-prevPos(:,1), intersectionPoint(2)-prevPos(:,2),prevPos.z]');
                            intersectionPoint(end) = obj.params.domainWidth/cos(theta);% z(prevPos.x,curPos.x,intersectionPoint(1),prevPos.y,curPos.y,intersectionPoint(2),prevPos.z,curPos.z);
                            domainNorm  = obj.GetDomainNormal(intersectionPoint);
                            cp      = [curPos.x curPos.y curPos.z];
                            pp      = [prevPos.x prevPos.y prevPos.z];
                            d       = (cp-pp)/sqrt(sum(cp-pp).^2);
                            r       = d-2*dot(d,domainNorm)*domainNorm;
                            n       = intersectionPoint+r*sqrt(sum(intersectionPoint-cp).^2)/sqrt(sum(r.^2));
                            curPos  = n;
                            prevPos = intersectionPoint;                            
                            obj.PlotReflection(pp,cp,intersectionPoint,n,domainNorm);
                        end
                        
                    elseif strcmpi(obj.params.domainShape,'twoPlates')
                        % The two plates are arranged vertically in parrallel
                        % the distance between the plates is
                        % obj.params.domainWidth
                        intersectionPoint = obj.FindIntersectionPoint(prevPos,curPos);% the intersection point is made sure to lay on the boundary
                        if ~isempty(intersectionPoint)
                            domainNorm        = obj.GetDomainNormal(intersectionPoint);
                            d                 = (curPos-prevPos)/sqrt(sum(curPos-prevPos).^2);
                            r                 = d-2*dot(d,domainNorm)*domainNorm;
                            n                 = intersectionPoint+r*sqrt(sum(intersectionPoint-curPos).^2)/sqrt(sum(r.^2));
                            obj.PlotReflection(prevPos,curPos,intersectionPoint,n,domainNorm);
                            curPos  = n;
                            prevPos = intersectionPoint;
                        end
                        
                    else
                        prevPos = pp;
                    end
                    
                    count = count+1;
                    if count>obj.params.maxReflectionsPerParticle
                        error('%s%s%s','Reflection count is bigger than',obj.params.maxReflectionsPerParticle,...
                            'Terminating');                        
                    end
                end
            end
        end
        
        function intersectionPoint = FindIntersectionPoint(obj,prevPos,curPos)

            if strcmpi(obj.params.domainShape,'sphere') || strcmpi(obj.params.domainShape,'cylinder')
                    % Find the intersection point
                    R    = obj.params.domainWidth;% domain radius
%                     L    = curPos-prevPos;        % direction of the particle
%                     d    = prevPos-obj.params.domainCenter;
%                     
%                     % the coefficient of the quadratic equation to solve
%                     A    = (sum(L.^2));
%                     B    = 2*dot(L,d);
%                     C    = (sum(d.^2))-R^2; 
%                     % the roots of the equation signify the length to go in the direction
%                     % vector, to reach the intersection point
%                     t(1) = (-B+sqrt(B^2-4*A*C))/(2*A);
%                     t(2) = (-B-sqrt(B^2-4*A*C))/(2*A);
                      %=== test ====
                     A  = prevPos;
                     B  = curPos; 
                     C  = (B-A);%./norm(B-A);
%                      dotCA = dot(C,A);
%                      dotCC = dot(C,C);
%                      dotAA = dot(A,A);
%                      t(1) = (-2*dotCA+sqrt(4*dotCA^2 -4*dotCC*(dotAA-R^2)))/(2*dotCC);
%                      t(2) = (-2*dotCA-sqrt(4*dotCA^2 -4*dotCC*(dotAA-R^2)))/(2*dotCC);
                     
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
                       % Make sure the intersection point is exactly on the
                       % boundary (roundoff and truncation error might lead
                       % to it being over the domain's boundary)
%                        er = (sqrt(sum(intersectionPoint.^2))-obj.params.domainWidth);
%                        if er~=0
%                            x = intersectionPoint(1);
%                            y = intersectionPoint(2);
%                            z = intersectionPoint(3);                           
%                            alpha(1) = (-2*(x+y+z)+sqrt(4*((x+y+z).^2) -12*((x.^2+y.^2+z.^2)-((R^2-2*er*R + er^2)))))/6;
%                            alpha(2) = (-2*(x+y+z)-sqrt(4*((x+y+z).^2) -12*((x.^2+y.^2+z.^2)-((R^2-2*er*R + er^2)))))/6;
%                            intersectionPoint = intersectionPoint+0.5*(min(real(alpha)));% zero out the last digit
% %                            intersectionPoint = (intersectionPoint/(sqrt(sum(intersectionPoint.^2)))) * obj.params.domainWidth;
%                        end
%                         er = (sqrt(sum(intersectionPoint.^2))-obj.params.domainWidth);
%                         if er~=0
%                             error('prob')
%                         end
                    else
                        intersectionPoint = [];
                    end
            elseif strcmpi(obj.params.domainShape,'twoPlates')
                    r  = obj.params.domainWidth;                                    
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
        
        function inIdx = InDomain(obj,vecIn)
            % Check if points are inside the domain 
            % the outut is a binary vector with 1 if the point is inside
            % the domain, 0 otherwise.
%             dimName       = {'x','y','z'};                                    
            numParticles = size(vecIn,1);
            inIdx        = true(numParticles,1);
            if strcmpi(obj.params.domainShape,'Sphere')
%                  v       = zeros(size(vecIn));
%             for dIdx = 1:obj.params.dimension
%                 v(:,dIdx)  = (vecIn.(dimName{dIdx}));                
%             end
%             v = vecIn
            % the vector norm
            n     = sqrt(sum(vecIn.^2,2));
            inIdx = n<=obj.params.domainWidth;
            
            elseif strcmpi(obj.params.domainShape,'cylinder')
%                  v       = zeros(numel(vecIn.x),obj.params.dimension);
%                 for dIdx = 1:2
%                     v(:,dIdx)  = (vecIn.(dimName{dIdx}));                
%                 end
                
            % the vector norm
            n     = sqrt(sum(vecIn.^2,2));
            inIdx = n<=obj.params.domainWidth;
            elseif strcmpi(obj.params.domainShape,'twoPlates')
                  inIdx = vecIn.x<obj.params.domainWidth;                 
            elseif strcmpi(obj.params.domainShape,'none')
                % Do nothing
                
            end
        end
        
        function domainNorm = GetDomainNormal(obj,intPoint)
           % Assuming the domain is a shpere 
            % the normals are facing inward
           if strcmpi(obj.params.domainShape,'sphere')
               domainNorm = obj.params.domainCenter-intPoint;
               domainNorm = domainNorm/sqrt(sum(domainNorm.^2));
               
           elseif strcmpi(obj.params.domainShape,'cylinder')
               domainNorm    = obj.params.domainCenter-intPoint;
               domainNorm(3) = 0;
               domainNorm = domainNorm/sqrt(sum(domainNorm.^2));
           elseif strcmpi(obj.params.domainShape,'twoPlates')
               if intPoint(1)< obj.params.domainWidth+1e-7 && intPoint(1)> obj.params.domainWidth-1e-7 
                   domainNorm = -[1,0 0];
               elseif intPoint(1)<-obj.params.domainWidth+1e-7 && intPoint(1)>-obj.params.domainWidth-1e-7
                   domainNorm = [1 0 0];
               end
           end
          
           
        end
                
        function PlotReflection(obj,prevPos,curPos,intPos,newPos,domainNorm)
            % the new pos is the point after reflection 
            
            % drw a line between the prev (in) point and cur (out) point
            line('XData',[prevPos(1) curPos(1)],...
                 'YData',[prevPos(2) curPos(2)],...
                 'ZData',[prevPos(3) curPos(3)],...
                 'Color','y',...
                 'Parent',obj.handles.graphical.mainAxes);
             
             % show intersection point
            line('XData',intPos(1),...
                 'YData',intPos(2),...
                 'ZData',intPos(3),...
                 'Marker','o',...
                 'MarkerEdgeColor','b',...
                 'lineStyle','none',...
                 'Color','m',...
                 'Parent',obj.handles.graphical.mainAxes);
             
             % circle the cur point
             line('XData',curPos(1),...
                 'YData',curPos(2),...
                 'ZData',curPos(3),...
                 'Marker','o',...
                 'MarkerEdgeColor','r',...
                 'lineStyle','none',...
                 'Parent',obj.handles.graphical.mainAxes);
             
             % circel prev point
             line('XData',prevPos(1),...
                 'YData',prevPos(2),...
                 'ZData',prevPos(3),...
                 'Marker','o',...
                 'MarkerEdgeColor','g',...
                 'lineStyle','none',...
                 'Parent',obj.handles.graphical.mainAxes);
             
             % draw a line between the reflection point and the new (in)
             % point
               line('XData',[intPos(1) newPos(1)],...
                 'YData',[intPos(2) newPos(2)],...
                 'ZData',[intPos(3) newPos(3)],...
                 'Marker','o',...                 
                 'MarkerEdgeColor','b',...
                 'lineStyle','-',...
                 'Color','c',...
                 'Parent',obj.handles.graphical.mainAxes);
             
             % show the domain normal
             hold on 
             quiver3(intPos(1),intPos(2),intPos(3),domainNorm(1), domainNorm(2), domainNorm(3),'Color','r')
%              quiver3(intPos(1),intPos(2),intPos(3),newPos(1), newPos(2), newPos(3),'Color','g')
             hold off
        end
    end
    
end