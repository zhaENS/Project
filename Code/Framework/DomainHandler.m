classdef DomainHandler<handle
    % At the moment the domain is parametrically defined 
    % TODO: allow insertion of a polygonal domain 
    properties 
        handles
        params
        points
    end
    
    events
    end
    
    methods 
        function obj = DomainHandler(domainParams)
            obj.params = domainParams;
            % class constructor 
%             obj.SetDefaultParams
%             
%             obj.SetInputParams(varargin{:})
            obj.handles.graphical.domain = [];
            obj.handles.graphical.light  = [];  
            obj.handles.graphical.mainAxes = obj.params.parentAxes;
            obj.ConstructDomain
                    
        end
        
        function SetDefaultParams(obj)
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
        
        function CreateControls(obj)
            obj.handles.graphical.mainFigure = figure('Tag','mainFigure',...
                'Units','normalized');
            
            obj.handles.graphical.mainAxes = axes('Parent',obj.handles.graphical.mainFigure,...
                'Tag','mainAxes',...
                'Units','normalized',...
                'NextPlot','ReplaceChildren');
            daspect([1 1 1])
        end
        
        function ConstructDomain(obj)
            if isempty(obj.handles.graphical.mainAxes)
                 obj.CreateControls
            else
%                 obj.handles.graphical.mainAxes = axesHandle;
            end
            
            if strcmpi(obj.params.domainShape,'sphere')
                [x,y,z] = sphere(20);
                
                obj.points.x = x*obj.params.domainWidth;
                obj.points.y = y*obj.params.domainWidth;
                obj.points.z = z*obj.params.domainWidth;
                
                if obj.params.dimension ==1
                    obj.points.y = zeros(size(y));
                    obj.points.z = zeros(size(z));
                elseif obj.params.dimension==2 
                    obj.points.z = zeros(size(z));
                end
            elseif strcmpi(obj.params.domainShape,'cylinder')
                % the cylinder z axis is pointing towatd [0 0 1];
                [x,y,z]      = cylinder(obj.params.domainWidth,20);
                obj.points.x = repmat(x,10,1);
                obj.points.y = repmat(y,10,1);
                obj.points.z = repmat(linspace(-obj.params.domainHeight/2,obj.params.domainHeight/2,...
                                    size(obj.points.x,1))',1,size(obj.points.x,2));  
                
               if obj.params.dimension ==1
                    obj.points.y = zeros(size(y));
                    obj.points.z = zeros(size(z));
                elseif obj.params.dimension==2 
                    obj.points.z = zeros(size(z));
                end
            elseif strcmpi(obj.params.domainShape,'twoPlates')
                   [obj.points.z,obj.points.y] = meshgrid(-obj.params.domainHeight:obj.params.domainHeight,...
                                                          -obj.params.domainHeight:obj.params.domainHeight,50);
                   obj.points.x    = obj.params.domainWidth*ones(size(obj.points.z));
                  obj.handles.graphical.mesh = mesh(obj.points.x,obj.points.y,obj.points.z,...
                                                     'FaceAlpha',0.3,...
                                                     'FaceColor',[0.7 0.9 0.2],...
                                                     'EdgeColor','none',...
                                                     'FaceLighting','phong',...
                                                     'Parent',obj.handles.graphical.mainAxes,...
                                                     'Tag','domain',...
                                                     'HandleVisibility','on',...
                                                     'Visible','off');
                 
                 [obj.points.z,obj.points.y] = meshgrid(-obj.params.domainHeight:obj.params.domainHeight,...
                                                         -obj.params.domainHeight:obj.params.domainHeight,50);
                  obj.points.x               = -obj.params.domainWidth*ones(size(obj.points.z));
                  
            elseif strcmpi(obj.params.domainShape,'none')
                    obj.points.x = [];
                    obj.points.y = [];
                    obj.points.z = [];

            end
            
            delete(obj.handles.graphical.domain);
            delete(obj.handles.graphical.light);
            
             m = mesh(obj.points.x,obj.points.y,obj.points.z,...
                     'FaceAlpha',0.3,...
                     'FaceColor',[0.7 0.9 0.2],...
                     'EdgeColor','none',...
                     'FaceLighting','phong',...
                     'Parent',obj.handles.graphical.mainAxes,...
                     'Tag','domain',...
                     'HandleVisibility','on',...
                     'Visible','off');
             
            l= light('Position',[1 1 1],...
                    'Style','local',...
                    'Tag','light',...
                    'Parent',obj.handles.graphical.mainAxes,...
                    'HandleVisibility','on',...
                    'Visible','off');  
             
            obj.handles.graphical.domain      = m;
            obj.handles.graphical.domainLight = l;
            
        end
        
        function ShowDomain(obj)
              f = findobj(get(obj.handles.graphical.mainAxes,'Children'),'Tag','domain');
              set(f,'Visible','on')
              l = findobj(get(obj.handles.graphical.mainAxes,'Children'),'Tag','light');
              set(l,'Visible','on');                        
        end
        
        function [prevPos,curPos,inFlag] = Reflect(obj,prevPos,curPos)%TODO: fix reflection for all domain shapes
%             disp('reflect')
            dimName  = {'x','y','z'};
            numBeads = numel(prevPos.x);
            inFlag   = true(numBeads,1);         
                count   = 1;
                while ~obj.InDomain(curPos)
                    % Reflect the particle
                    if strcmpi(obj.params.domainShape,'none')
                        % Do nothing
                    elseif strcmpi(obj.params.domainShape,'sphere')
                        intersectionPoint = obj.FindIntersectionPoint([prevPos.x prevPos.y prevPos.z],[curPos.x curPos.y curPos.z]);
                        if obj.InDomain(prevPos) && ~obj.InDomain(curPos) && isempty(intersectionPoint) %%|| sqrt(sum(intersectionPoint.^2))~=obj.params.domainWidth)
                            error('something is wrong with the intersection point function')
                        end
                        if ~isempty(intersectionPoint)
                            domainNorm  = obj.GetDomainNormal(intersectionPoint);
                            
                            % Calculate reflected ray
                            pp  = [prevPos.x prevPos.y prevPos.z];
                            cp  = [curPos.x curPos.y curPos.z];
                            di  = (intersectionPoint-pp)/sqrt(sum(intersectionPoint-pp).^2);
                            ds  = -(2*dot(domainNorm,di)*domainNorm-di);
                            t   = sqrt(sum(cp-intersectionPoint).^2);
                            n   = intersectionPoint+ t*ds;% the new position
                            
%                             obj.PlotReflection(pp,cp,intersectionPoint,n,domainNorm);
                            biasNorm = (n-intersectionPoint)/norm(n-intersectionPoint);
                            %
                            %
                            for dIdx = 1:obj.params.dimension
                                curPos.(dimName{dIdx})  = n(dIdx);
                                % To avoid numerical error, move the prev point
                                % slightly on the vector between the new point and the
                                % intersectionPoint
                                prevPos.(dimName{dIdx}) = intersectionPoint(dIdx)+(1e-15)*biasNorm(dIdx);
                                
                            end
                        else
                            disp('no intersection point')
                        end
                        
                    elseif strcmpi(obj.params.domainShape,'cylinder');
                        % find the intersection in 2D first
                        intersectionPoint = obj.FindIntersectionPoint([prevPos.x prevPos.y 0],[curPos.x curPos.y 0]);
                        
                        if ~isempty(intersectionPoint)
                            theta =subspace([curPos.x-prevPos.x curPos.y-prevPos.y curPos.z-prevPos.z]',[intersectionPoint(1)-prevPos.x, intersectionPoint(2)-prevPos.y,prevPos.z]');
                            intersectionPoint(end) = obj.params.domainWidth/cos(theta);% z(prevPos.x,curPos.x,intersectionPoint(1),prevPos.y,curPos.y,intersectionPoint(2),prevPos.z,curPos.z);
                            domainNorm  = obj.GetDomainNormal(intersectionPoint);
                            cp      = [curPos.x curPos.y curPos.z];
                            pp      = [prevPos.x prevPos.y prevPos.z];
                            d       = (cp-pp)/sqrt(sum(cp-pp).^2);
                            r       = d-2*dot(d,domainNorm)*domainNorm;
                            n       = intersectionPoint+r*sqrt(sum(intersectionPoint-cp).^2)/sqrt(sum(r.^2));
                            curPos  = struct('x',n(1),'y',n(2),'z',n(3));
                            prevPos = struct('x',intersectionPoint(1),'y',intersectionPoint(2),'z',intersectionPoint(3));
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
                            newPos  = struct('x',[],'y',[],'z',[]);
                            prevPos = struct('x',[],'y',[],'z',[]);
                            for dIdx = 1:obj.params.dimension
                                newPos.(dimName{dIdx})  = n(dIdx);
                                prevPos.(dimName{dIdx}) = intersectionPoint(dIdx);
                            end
                        end
                        
                    else
                        %                 newPos  = cp;
                        prevPos = pp;
                    end
                    count = count+1;
                    if count>10
                        disp('reflection count is bigger than 10. Terminating');
                        error('too many reflection iterations')
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
                     dotCA = dot(C,A);
                     dotCC = dot(C,C);
                     dotAA = dot(A,A);
                     t(1) = (-2*dotCA+sqrt(4*dotCA^2 -4*dotCC*(dotAA-R^2)))/(2*dotCC);
                     t(2) = (-2*dotCA-sqrt(4*dotCA^2 -4*dotCC*(dotAA-R^2)))/(2*dotCC);
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
            dimName = {'x','y','z'};                                    
           
            inIdx = ones(numel(vecIn.x),1);
            if strcmpi(obj.params.domainShape,'Sphere')
                 v       = zeros(numel(vecIn.x),obj.params.dimension);
            for dIdx = 1:obj.params.dimension
                v(:,dIdx)  = (vecIn.(dimName{dIdx}));                
            end
            % the vector norm
            n     = sqrt(sum(v.^2,2));
            inIdx = n<=obj.params.domainWidth;
            
            elseif strcmpi(obj.params.domainShape,'cylinder')
                 v       = zeros(numel(vecIn.x),obj.params.dimension);
                for dIdx = 1:2
                    v(:,dIdx)  = (vecIn.(dimName{dIdx}));                
                end
            % the vector norm
            n     = sqrt(sum(v.^2,2));
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