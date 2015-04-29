classdef SimulationFrameworkGraphics<handle
    % This class controls the graphical features of the simulation
    % framework
    % TODO: graphics should operate offline after simulations are done
    properties
        runSimulation
        handles
        mainAxes        
        chainColors
        params % holds the framework parameters
%         chain = struct('bead',[],'connectors',[])
    end
    
    methods
        function obj = SimulationFrameworkGraphics(frameworkHandle)
            obj.handles.framework = frameworkHandle;
            obj.params            = frameworkHandle.params;           
        end
        
        function CreateControls(obj)
           if obj.params.simulator.showSimulation
             % Main figure          
            obj.handles.graphical.mainFigure = figure('Units','normalized',...
                                                      'Tag','mainFigure',...
                                                      'ToolBar','none',...
                                                      'MenuBar','none',...                                                     
                                                      'BusyAction','queue');
                                                      
            % Main axes                                      
            obj.handles.graphical.mainAxes  = axes('Units','normalized',...
                                                   'Parent',obj.handles.graphical.mainFigure,...
                                                   'NextPlot','ReplaceChildren',...
                                                   'Color',[0.2 0.1 0.23],...
                                                   'Box','on',...
                                                   'XTick',[],...
                                                   'YTick',[],...
                                                   'ZTick',[]);      
         
           
           obj.handles.graphical.stop   = uicontrol('Units','normalized',...
               'Parent',obj.handles.graphical.mainFigure,...
               'String','Stop',...
               'Tag','StartStopButton',...
               'Position',[0.05 0.05, 0.1, 0.05],...
               'Callback',@obj.StartStop);
           
           obj.handles.graphical.stop   = uicontrol('Units','normalized',...
               'Parent',obj.handles.graphical.mainFigure,...
               'String','Next',...
               'Tag','nextButton',...
               'Position',[0.15 0.05, 0.1, 0.05],...
               'Callback',@obj.Step);
                           
                % show domain  
                obj.ConstructDomain;                
                
                % show chains 
                obj.InitializeChainGraphics
                
                cameratoolbar
            end
        end
        
        function ConstructDomain(obj)
            
            if strcmpi(obj.params.domain.domainShape,'sphere')
                [x,y,z] = sphere(20);
                
                obj.handles.graphical.domain.points.x = x*obj.params.domain.domainWidth;
                obj.handles.graphical.domain.points.y = y*obj.params.domain.domainWidth;
                obj.handles.graphical.domain.points.z = z*obj.params.domain.domainWidth;
                
                if obj.params.domain.dimension ==1
                    obj.handles.graphical.domain.points.y = zeros(size(y));
                    obj.handles.graphical.domain.points.z = zeros(size(z));
                elseif obj.params.domain.dimension==2 
                    obj.handles.graphical.domain.points.z = zeros(size(z));
                end
            elseif strcmpi(obj.params.domain.domainShape,'cylinder')
                % the cylinder z axis is pointing towatd [0 0 1];
                [x,y,z]      = cylinder(obj.params.domain.domainWidth,20);
                obj.handles.graphical.domain.points.x = repmat(x,10,1);
                obj.handles.graphical.domain.points.y = repmat(y,10,1);
                obj.handles.graphical.domain.points.z = repmat(linspace(-obj.params.domain.domainHeight/2,obj.params.domain.domainHeight/2,...
                                    size(obj.handles.graphical.domain.points.x,1))',1,size(obj.handles.graphical.domain.points.x,2));  
                
               if obj.params.domain.dimension ==1
                    obj.handles.graphical.domain.points.y = zeros(size(y));
                    obj.handles.graphical.domain.points.z = zeros(size(z));
                elseif obj.params.dimension==2 
                    obj.handles.graphical.domain.points.z = zeros(size(z));
                end
            elseif strcmpi(obj.params.domain.domainShape,'twoPlates')
                   [obj.handles.graphical.domain.points.z,obj.handles.graphical.domain.points.y] = ...
                          meshgrid(-obj.params.domain.domainHeight:obj.params.domainHeight,...
                           -obj.params.domain.domainHeight:obj.params.domain.domainHeight,50);
                       
                   obj.handles.graphical.domain.points.x    = obj.params.domain.domainWidth*ones(size(obj.handles.graphical.domain.points.z));
                   obj.handles.graphical.domain.mesh = mesh(obj.handles.graphical.domain.points.x,...
                                                            obj.handles.graphical.domain.points.y,...
                                                            obj.handles.graphical.domain.points.z,...
                                                             'FaceAlpha',0.3,...
                                                             'FaceColor',[0.7 0.9 0.2],...
                                                             'EdgeColor','none',...
                                                             'FaceLighting','phong',...
                                                             'Parent',obj.handles.graphical.mainAxes,...
                                                             'Tag','domain',...
                                                             'HandleVisibility','on',...
                                                             'Visible','on');
                 
                 [obj.handles.graphical.domain.points.z,obj.handles.graphical.domain.points.y] = ...
                                             meshgrid(-obj.params.domainHeight:obj.params.domainHeight,...
                                                -obj.params.domain.domainHeight:obj.params.domainHeight,50);
                  obj.handles.graphical.domain.points.x  = -obj.params.domain.domainWidth*ones(size(obj.handles.graphical.domain.points.z));
                  
            elseif strcmpi(obj.params.domain.domainShape,'open')
                    obj.handles.graphical.domain.points.x = [];
                    obj.handles.graphical.domain.points.y = [];
                    obj.handles.graphical.domain.points.z = [];
            else 
                error('unsupported domain type');
            end            
            
             m = mesh(obj.handles.graphical.domain.points.x,...
                     obj.handles.graphical.domain.points.y,...
                     obj.handles.graphical.domain.points.z,...
                     'FaceAlpha',0.3,...
                     'FaceColor',[0.7 0.9 0.2],...
                     'EdgeColor','none',...
                     'FaceLighting','phong',...
                     'Parent',obj.handles.graphical.mainAxes,...
                     'Tag','domain',...
                     'HandleVisibility','on',...
                     'Visible','on');
             
            l= light('Position',[1 1 1],...
                    'Style','local',...
                    'Tag','light',...
                    'Parent',obj.handles.graphical.mainAxes,...
                    'HandleVisibility','on',...
                    'Visible','on');  
             
            obj.handles.graphical.domain      = m;
            obj.handles.graphical.domainLight = l;
            
      end       
        
        function ShowSimulation(obj)
         if obj.params.simulator.showSimulation 
            obj.UpdateChainPosition
         end
        end
                
        function UpdateChainPosition(obj)
               %TODO: plotting should be able to take into account
               %composite objects and not only chains 

                numObjects   = obj.handles.framework.objectManager.numObjects;
                curPos       = obj.handles.framework.objectManager.curPos;% all positions
                
                for oIdx = 1:numObjects;

                    % get the members of the object 
                    numMembers = obj.handles.framework.objectManager.map.GetObjectCount(oIdx);

                    member = obj.handles.framework.objectManager.map.object(oIdx).members; % member index in the general list

                    color =obj.chainColors(member(1),:);% transform to the first member's color 
                    for mIdx = 1:numMembers                        
                        inds = obj.handles.framework.objectManager.map.GetMemberInds(oIdx,mIdx);
                        c    = curPos(inds,:);
                        
                    % Plot beads                   
                    set(obj.handles.graphical.chain(member(mIdx)).beads,...
                        'XData',c(:,1),...
                        'YData',c(:,2),...
                        'ZData',c(:,3),...
                        'MarkerEdgeColor',color,...
                        'MarkerFaceColor',color,...
                        'Color',color);
%                         'Visible','on');
                   
%                     Show connectors between beads 
%                         par = chainParams{oIdx};
%                         cm  = [par(:).connectedBeads];
%                         cm  = obj.handles.framework.objectManager.GetConnectedParticles(oIdx);            
% %                      if ~isempty(cm)
%                         set(obj.handles.graphical.chain(member(mIdx)).connectors,...
%                                 'XData',[c(cm(:,1),1);c(cm(:,2),1)],...
%                                 'YData',[c(cm(:,1),2);c(cm(:,2),2)],...
%                                 'ZData',[c(cm(:,1),3);c(cm(:,2),3)],...
%                                 'Color','y',...
%                                 'LineWidth',4,...
%                                 'Tag','connector');
                                
                     end
%                      else
%                        % Create a phantom connector, a place holder for future
%                        % connections                     
%                      set(obj.handles.graphical.chain(oIdx).connectors,...
%                                'XData',[0;0],...
%                                'YData',[0;0],...
%                                'ZData',[0;0],...
%                                'Tag','connector',...
%                                'Color','w',...
%                                'Visible','off');
                  
%                      end   
%                     end
%                 end
                end
%             hide all non-active objects
%             hideObj = setdiff(1:size(obj.chainColors,1),1:numObjects);
%             for hIdx = 1:numel(hideObj)
% %                 member = obj.handles.framework.objectManager.map.object(hideObj(hIdx)).members 
% %                 for mIdx = 1:numel(member)
%                  set(obj.handles.graphical.chain(hideObj(hIdx)).connectors,...
%                                'XData',[0;0],...
%                                'YData',[0;0],...
%                                'ZData',[0;0],...
%                                'Tag','connector',...
%                                'Color','w',...
%                                'Visible','off');
%                 set(obj.handles.graphical.chain(hideObj(hIdx)).beads,...
%                         'XData',0,...
%                         'YData',0,...
%                         'ZData',0,...
%                         'Visible','off');
% %                 end
%             end
            drawnow
        end
        
        function Step(obj,varargin)
            obj.handles.framework.Step
        end
        
        function StartStop(obj,buttonHandle,varargin)
            if obj.handles.framework.runSimulation 
                set(buttonHandle,'String','Continue')
                obj.handles.framework.runSimulation = false;
            else 
                set(buttonHandle,'String','Stop')
                obj.handles.framework.runSimulation = true;
                obj.handles.framework.Run
            end
        end
        
        function InitializeChainGraphics(obj)%TODO: discard the matrix representation for connectivity
            % set graphics
           numObj           = obj.handles.framework.objectManager.numObjects;
           [prevPos,curPos] = obj.handles.framework.objectManager.GetPosition(1:numObj);
           objParams        = obj.handles.framework.objectManager.GetObjectParameters(1:numObj);
           obj.handles.graphical.chain= struct('beads',[],'connectors',[]);
           
           % choose chain colors 
           l = linspace(0.4,0.9,obj.params.simulator.numChains);
          
           for oIdx = 1:numObj
               r = zeros(1,3);
               for rIdx = 1:3
                k       = randperm(obj.params.simulator.numChains);
                r(rIdx) = k(1);
               end
               obj.chainColors(oIdx,1:3) = [l(r(1)),l(r(2)),l(r(3))];
           end
           
           
           
        % Create chains 
            for oIdx = 1:numObj
                
                if ~(objParams{oIdx}.springForce) && ~objParams{oIdx}.bendingElasticityForce
                    lineStyle = 'none';
                else
                    lineStyle = '-';
                end
                
                % Plot beads
                p     = prevPos{oIdx};
                c     = curPos{oIdx};
                obj.handles.graphical.chain(oIdx).beads = line(...
                    'XData',p(:,1),...                                   
                    'YData',p(:,2),...
                    'ZData',p(:,3),...
                    'Marker','o',...
                    'LineStyle',lineStyle,...
                    'Color','w',...
                    'Color',obj.chainColors(oIdx,:),...
                    'Tag','chain',...
                    'MarkerSize',6,...
                    'MarkerFaceColor',obj.chainColors(oIdx,:),...
                    'Parent',obj.handles.graphical.mainAxes,...
                    'Visible','on');
                
                % Plot connectors                  
                cm = objParams{oIdx}.connectedBeads;                               
                   
                if ~isempty(cm)
                   obj.handles.graphical.chain(oIdx).connectors=line(...
                               [c(cm(:,1),1);c(cm(:,2),1)],...
                               [c(cm(:,1),2);c(cm(:,2),2)],...
                               [c(cm(:,1),3);c(cm(:,2),3)],...
                                'Tag','connector',...
                                'Color','y',...
                                'LineWidth',4,...
                                'Visible','on');
                else 
                  % Create a phantom connector, a place holder for future
                % connections 
                   obj.handles.graphical.chain(oIdx).connectors=line(...
                                                                   [0;0],...
                                                                   [0;0],...
                                                                   [0;0],...
                                                                   'Tag','connector',...
                                                                   'Color','w',...
                                                                   'Visible','off');  
                end              
            end 
        end
        
        function SetAxesLimits(obj)% fix
            if strcmpi(obj.params.domainShape,'open')
            bias = obj.params.domainWidth+2;
            set(obj.handles.graphical.mainAxes,'XLim',[obj.chainsCenterOfMass.x-bias obj.chainsCenterOfMass.x+bias],...
                'YLim',[obj.chainsCenterOfMass.y-bias obj.chainsCenterOfMass.y+bias],...
                'ZLim',[obj.chainsCenterOfMass.z-bias obj.chainsCenterOfMass.z+bias]);
            else
            bias = obj.params.domainWidth+5;
            
            set(obj.handles.graphical.mainAxes,...
                'XLim',[0-bias, 0+bias],...
                'YLim',[0-bias, 0+bias],...
                'ZLim',[0-bias, 0+bias]);
            end
             daspect([1 1 1])
%              drawnow
        end
        
    end
end