classdef ChainDynamicsPlayer<handle
    properties
        params
        handles
        chain
        affineBeads
        connectedBeads        
    end
    
    properties (Access = private)
        play = false;
        prevStep  = 1;
    end
    
    methods
        function obj = ChainDynamicsPlayer(rouseChain)
            % construct the class with rouseChain (from SimpleRouse Class)
            if exist('rouseChain','var')
                obj.handles.classes.rouseChain = rouseChain;
                obj.affineBeads      = rouseChain.params.affineBeadsNum;
                obj.connectedBeads   = rouseChain.params.connectedBeads; 
                obj.chain            = rouseChain.savedPosition;
                obj.params.dimension = rouseChain.params.dimension;
                obj.params.numSteps  = rouseChain.step;
                obj.params.numBeads  = rouseChain.params.numBeads;
                obj.params.connectorRadius = 4;
                obj.params.connectorColor   = [0.1,0.8,0.1];
                obj.params.beadRadius = 0.1;
                obj.params.beadColor  = [0.1 0.2, 0.7];
                obj.params.affineBeadColor  = [0.1, 0.5 0.6];
                
                obj.CreateControls
                obj.Display(1);% display first step 
            end
        end
        
        function CreateControls(obj)
            %
            obj.handles.graphical.mainFigure = figure('Units','norm');
            % create main axes
            obj.handles.graphical.mainAxes  = axes('Units','norm',...
                                                    'NextPlot','Add',...
                                                    'Parent',obj.handles.graphical.mainFigure,...
                                                    'XTick',[],...
                                                    'YTick',[],...
                                                    'ZTick',[]);
            box on 
            
            % add light object 
            obj.handles.graphical.light = light('Parent',obj.handles.graphical.mainAxes);
            
            % add camera tool bar
            obj.handles.graphical.cameraToolBar = cameratoolbar(obj.handles.graphical.mainFigure);
            
            % set axes limits 
            obj.SetAxesLimits
            
            % construct step slider
            obj.handles.graphical.stepSlider = uicontrol('Parent',obj.handles.graphical.mainFigure,...
                'Style','slider',...
                'Max',obj.params.numSteps,...
                'Min', 1,...
                'Units','norm',...
                'Value',1,...
                'Position',[0.95, 0, 0.05, 0.5],...
                'SliderStep',[1/(obj.params.numSteps-1), 1/(obj.params.numSteps-1)],...
                'Callback',@obj.MoveSlider);
            
            % construct play button 
            obj.handles.graphical.playButton = uicontrol('Parent',obj.handles.graphical.mainFigure,...
                'Units','norm',...
                'Position',[0.95, 0.55,0.05,0.05],...
                'String','Play',...
                'Callback',@obj.Play);
            
%             % initial chain display
                         
            % add the beads as spheres
            for bIdx = 1:obj.params.numBeads
                [x,y,z] = sphere(10);
                x = obj.params.beadRadius*x+obj.chain(bIdx,1,1);
                y = obj.params.beadRadius*y+obj.chain(bIdx,2,1);
                z = obj.params.beadRadius*z+obj.chain(bIdx,3,1);
                if ismember(bIdx,obj.affineBeads(:))
                    beadColor =  obj.params.affineBeadColor ;
                else
                    beadColor =  obj.params.beadColor ;
                end
                
                obj.handles.graphical.beads(bIdx) = mesh(x,y,z,...
                                                        'EdgeColor','none',...
                                                        'Facecolor',beadColor,...
                                                        'FaceLighting','phong',...
                                                        'Parent',obj.handles.graphical.mainAxes);
            end
            
            
            
            obj.handles.graphical.chain = line('XData',obj.chain(:,1,1),...
                'YData',obj.chain(:,2,1),...
                'ZData',obj.chain(:,3,1),...
                'Parent',obj.handles.graphical.mainAxes,...
                'LineWidth',obj.params.connectorRadius,...
                'Color',obj.params.connectorColor,...
                'Marker','none');
            
%             for aIdx = 1:size(obj.affineBeads,1)
%                 obj.handles.graphical.affineBeads(aIdx) = line('XData',obj.chain(obj.affineBeads(aIdx,:),1,1),...
%                     'YData',obj.chain(obj.affineBeads(aIdx,:),2,1),...
%                     'ZData',obj.chain(obj.affineBeads(aIdx,:),3,1),...
%                     'Parent',obj.handles.graphical.mainAxes,...
%                     'Marker','o',...
%                     'LineStyle','none',...
%                     'MarkerFaceColor',rand(1,3));
%             end
            
            for cIdx = 1:size(obj.connectedBeads,1)
                obj.handles.graphicsl.connectedBeads(cIdx) = line(...
                    'XData',[obj.chain(obj.connectedBeads(cIdx,1),1,1) obj.chain(obj.connectedBeads(cIdx,2),1,1)],...
                    'YData',[obj.chain(obj.connectedBeads(cIdx,1),2,1),obj.chain(obj.connectedBeads(cIdx,2),2,1)],...
                    'ZData',[obj.chain(obj.connectedBeads(cIdx,1),3,1),obj.chain(obj.connectedBeads(cIdx,2),3,1)],...
                    'Parent',obj.handles.graphical.mainAxes,...
                    'Marker','none',...
                    'LineStyle','-',...
                    'Color',obj.params.connectorColor,...
                    'LineWidth',obj.params.connectorRadius);
            end
            
        end
        
        function Display(obj,stepNum)
            
            % Move chain connectors
            set(obj.handles.graphical.chain,'XData',obj.chain(:,1,stepNum),...
                'YData',obj.chain(:,2,stepNum),...
                'ZData',obj.chain(:,3,stepNum),...
                'Parent',obj.handles.graphical.mainAxes);
            
%             for aIdx = 1:size(obj.affineBeads,1)
%                 set(obj.handles.graphical.affineBeads(aIdx),...
%                     'XData',obj.chain(obj.affineBeads(aIdx,:),1,stepNum),...
%                     'YData',obj.chain(obj.affineBeads(aIdx,:),2,stepNum),...
%                     'ZData',obj.chain(obj.affineBeads(aIdx,:),3,stepNum),...
%                     'Parent',obj.handles.graphical.mainAxes,...
%                     'Marker','o',...
%                     'LineStyle','none');
%             end
            
            % move beads spheres
             prevX = obj.chain(:,1,obj.prevStep);
             prevY = obj.chain(:,2,obj.prevStep);
             prevZ = obj.chain(:,3,obj.prevStep);                        
            for bIdx = 1:obj.params.numBeads
                set(obj.handles.graphical.beads(bIdx),'XData',get(obj.handles.graphical.beads(bIdx),'XData')-prevX(bIdx)+obj.chain(bIdx,1,stepNum),...
                    'YData',get(obj.handles.graphical.beads(bIdx),'YData')-prevY(bIdx)+obj.chain(bIdx,2,stepNum),...
                    'ZData',get(obj.handles.graphical.beads(bIdx),'ZData')-prevZ(bIdx)+obj.chain(bIdx,3,stepNum));
            end
            
            
            for cIdx = 1:size(obj.connectedBeads,1)
                 set(obj.handles.graphicsl.connectedBeads(cIdx),....
                    'XData',[obj.chain(obj.connectedBeads(cIdx,1),1,stepNum) obj.chain(obj.connectedBeads(cIdx,2),1,stepNum)],...
                    'YData',[obj.chain(obj.connectedBeads(cIdx,1),2,stepNum),obj.chain(obj.connectedBeads(cIdx,2),2,stepNum)],...
                    'ZData',[obj.chain(obj.connectedBeads(cIdx,1),3,stepNum),obj.chain(obj.connectedBeads(cIdx,2),3,stepNum)],...
                    'Parent',obj.handles.graphical.mainAxes,...
                    'Marker','none',...
                    'LineStyle','-',...
                    'LineWidth',obj.params.connectorRadius,...
                    'MarkerFaceColor','b');
            end
            drawnow
        end
        
        function Play(obj,bHandle,varargin)
            cStep = round(get(obj.handles.graphical.stepSlider,'Value'));
            if strcmpi(get(bHandle,'String'),'Play')
                set(bHandle,'String','Pause')
                obj.play = true;
                while obj.play
                    
                    obj.prevStep = cStep;
                    cStep = cStep+1;
                    set(obj.handles.graphical.stepSlider,'Value',cStep)
                     
                    obj.Display(cStep)
                    if cStep>=obj.params.numSteps
                        obj.play = false;
                    end
                end
            else
                obj.play = false;
                set(bHandle,'String','Play')
            end
        end
        
        function SetAxesLimits(obj)
            x =obj.chain(:,1,:);
            x = x(:);
            y = obj.chain(:,2,:);
            y = y(:);
            z = obj.chain(:,3,:);
            z = z(:);
            maxX = max(x);
            minX = min(x);
            maxY = max(y);
            minY = min(y);
            maxZ = max(z);
            minZ = min(z);
            set(obj.handles.graphical.mainAxes,'XLim',[minX,maxX],...
                'YLim',[minY, maxY],...
                'ZLim',[minZ,maxZ]);                        
        end
        
        function MoveSlider(obj,sliderHandle,varargin)
            sv  = ceil(get(sliderHandle,'Value'));
            obj.Display(sv);
            obj.prevStep = sv;
        end
    end
    
end