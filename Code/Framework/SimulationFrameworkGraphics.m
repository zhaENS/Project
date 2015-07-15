classdef SimulationFrameworkGraphics<handle
    % This class controls the graphical features of the simulation
    % framework
    % TODO: graphics should operate offline after simulations are done
    % TODO: include domain center in domain construction
    % TODO: incorporate Histone graphics 
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
                                
                obj.handles.graphical.startStopButton   = uicontrol('Units','normalized',...
                    'Parent',obj.handles.graphical.mainFigure,...
                    'String','Stop',...
                    'Tag','StartStopButton',...
                    'Position',[0.05 0.05, 0.1, 0.05],...
                    'Callback',@obj.StartStop);
                
                obj.handles.graphical.nextButton   = uicontrol('Units','normalized',...
                    'Parent',obj.handles.graphical.mainFigure,...
                    'String','Next',...
                    'Tag','nextButton',...
                    'Position',[0.15 0.05, 0.1, 0.05],...
                    'Callback',@obj.Step);
                
                obj.handles.graphical.snapshotButton = uicontrol('Units','normalized',...
                    'Parent',obj.handles.graphical.mainFigure,...
                    'String','snapshot',...
                    'Tag','snapshotButton',...
                    'Position',[0.05, 0.15,0.15,0.05],...
                    'Callback',@obj.Snapshot);
                
                % show domain
                obj.ConstructDomain;
                
                % show chains
                obj.InitializeChainGraphics
                
                cameratoolbar
            end
        end
        
        function ConstructDomain(obj)
            numDomains = numel(obj.params.domain);
            set(obj.handles.graphical.mainAxes,'NextPlot','Add')
            for dIdx = 1:numDomains
                dp = obj.params.domain(dIdx);
                if strcmpi(dp.domainShape,'sphere')
                    [x,y,z] = sphere(20);
                    
                    points.x = x*dp.domainWidth+ dp.domainCenter(1);
                    points.y = y*dp.domainWidth+ dp.domainCenter(2);
                    points.z = z*dp.domainWidth+ dp.domainCenter(3);
                    
                    if dp.dimension ==1
                        points.y = zeros(size(y));
                        points.z = zeros(size(z));
                    elseif dp.dimension==2
                        points.z = zeros(size(z));
                    end
                elseif strcmpi(dp.domainShape,'cylinder')
                    % the cylinder z axis is pointing towatd [0 0 1];
                    [x,y,z]      = cylinder(dp.domainWidth,20);
                    points.x = repmat(x,10,1);
                    points.y = repmat(y,10,1);
                    points.z = repmat(linspace(-dp.domainHeight/2,dp.domainHeight/2,...
                        size(points.x,1))',1,size(points.x,2));
                    
                    if dp.dimension ==1
                        points.y = zeros(size(y));
                        points.z = zeros(size(z));
                    elseif dp.dimension==2
                        points.z = zeros(size(z));
                    end
                elseif strcmpi(dp.domainShape,'twoPlates')
                    [points.z,points.y] = ...
                        meshgrid(-dp.domainHeight:dp.domainHeight,...
                        -dp.domainHeight:dp.domainHeight,50);
                    
                    points.x    = dp.domainWidth*ones(size(points.z));
                    
                    [points.z,points.y] = ...
                        meshgrid(-dp.domainHeight:dp.domainHeight,...
                        -dp.domainHeight:dp.domainHeight,50);
                    points.x  = -dp.domainWidth*ones(size(points.z));
                    
                elseif strcmpi(dp.domainShape,'open')
                    points.x = [];
                    points.y = [];
                    points.z = [];
                else
                    error('unsupported domain type');
                end
                
                % assign the points handles
                obj.handles.graphical.domain(dIdx).points = points;
                
                % create a mesh
                obj.handles.graphical.domain(dIdx).mesh = mesh(obj.handles.graphical.mainAxes,points.x,points.y,points.z,...
                                'FaceAlpha',0.3,...
                                'FaceColor',[0.7 0.9 0.2],...
                                'EdgeColor','none',...
                                'FaceLighting','phong',...                                
                                'Tag','domain',...
                                'HandleVisibility','on',...
                                'Visible','on');
                

            end
            
            l= light('Position',[1 1 1],...
                'Style','infinite',...
                'Tag','light',...
                'Parent',obj.handles.graphical.mainAxes,...
                'HandleVisibility','on',...
                'Visible','on');
            
            
            set( obj.handles.graphical.domain(dIdx).mesh,'FaceColor',[0.7 0.2, 0.4],'FaceAlpha',0.3);% temp
            
            obj.handles.graphical.domainLight = l;
            set(obj.handles.graphical.mainAxes,'NextPlot','ReplaceChildren')
            daspect(obj.handles.graphical.mainAxes,[1 1 1])
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
                
                member = obj.handles.framework.objectManager.map.GetObjectMembers(oIdx); % member index in the general list
                
                color = obj.chainColors(member(1),:);% transform to the first member's color
                tubeFlag = false;
                for mIdx = 1:numMembers
                    inds = obj.handles.framework.objectManager.map.GetMemberInds(oIdx,mIdx);
                    c    = curPos(inds,:);
                    if tubeFlag
                        hold on
                        [x,y,z,sh] = tubeplot([c(:,1),c(:,2),c(:,3)]',0.1);
                        
                        %                                 sh = surf(obj.handles.graphical.mainAxes,x,y,z);
                        obj.handles.graphical.chain(member(mIdx)).beads = sh;
                        hold off
                        drawnow
                    else
                        
                        % Plot beads
                        set(obj.handles.graphical.chain(member(mIdx)).beads,...
                            'XData',c(:,1),...
                            'YData',c(:,2),...
                            'ZData',c(:,3),...
                            'MarkerEdgeColor',color,...
                            'MarkerFaceColor',color,...
                            'Color',color);
                    end
                    %                     Show connectors between beads
                    %                     cm = obj.handles.framework.objectManager.GetObjectConnectedParticles(oIdx,'offDiagonals');
                    %
                    %                      if ~isempty(cm)
                    %                         set(obj.handles.graphical.chain(member(mIdx)).connectors,...
                    %                                 'XData',[curPos(cm(:,1),1);curPos(cm(:,2),1)],...
                    %                                 'YData',[curPos(cm(:,1),2);curPos(cm(:,2),2)],...
                    %                                 'ZData',[curPos(cm(:,1),3);curPos(cm(:,2),3)],...
                    %                                 'Color','y',...
                    %                                 'LineWidth',5,...
                    %                                 'Tag','connector',...
                    %                                 'Visible','on');
                    %                      else
                    % Create a phantom connector, a place holder for future
                    %                        % connections
                    %                      set(obj.handles.graphical.chain(member(mIdx)).connectors,...
                    %                                'XData',[0;0],...
                    %                                'YData',[0;0],...
                    %                                'ZData',[0;0],...
                    %                                'Tag','connector',...
                    %                                'Color','w',...
                    %                                'Visible','off');
                    
                    %                      end
                end
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
        
        function Snapshot(obj,varargin)
            
            c       = clock;
            % assign a name according to the clock 
            uName   = ['snapshot_' num2str(c(2:end))];
            [~,~,~] = mkdir(obj.params.plotHandler.galleryFolderPath);
            fr      = obj.handles.framework;
            
            % hide control for the saving 
            set([obj.handles.graphical.startStopButton,...
                 obj.handles.graphical.snapshotButton,...
                 obj.handles.graphical.nextButton],'Visible','off');
            
             t = title(obj.handles.graphical.mainAxes,...
                      ['Step: ',num2str(fr.simulationData(fr.batchRound,fr.simulationRound).step),...
                       ' Time: ',sprintf('%f',fr.simulationData(fr.batchRound,fr.simulationRound).step*fr.params.simulator.dt)]);
                                              
            saveas(obj.handles.graphical.mainFigure,fullfile(obj.params.plotHandler.galleryFolderPath,uName),'png');
            hgsave(obj.handles.graphical.mainFigure,fullfile(obj.params.plotHandler.galleryFolderPath,[uName,'.fig']),'-v7.3');
             hgsave(obj.handles.graphical.mainFigure,fullfile(obj.params.plotHandler.galleryFolderPath,[uName,'.eps']),'-v7.3');
            set([obj.handles.graphical.startStopButton,...
                 obj.handles.graphical.snapshotButton,...
                 obj.handles.graphical.nextButton],'Visible','on')
           % delete(t);
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
                tubeFlag = false;
                if tubeFlag
                    hold on
                    [~,~,~,sh] = tubeplot([c(:,1),c(:,2),c(:,3)]',0.1);
                    obj.handles.graphical.chain(oIdx).beads = sh;
                    %                     surf(obj.handles.graphical.mainAxes,get(sh,'XData'),get(sh,'YData'),get(sh,'ZData'));
                    hold off
                else
                    
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
                end
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