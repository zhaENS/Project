classdef SimulationFrameworkGraphics<handle
    properties        
        runSimulation
        handles
        mainAxes        
        chainColors
        params % holds the framework parameters
        chain = struct('bead',[],'connectors',[])
    end
    
    methods
        function obj = SimulationFrameworkGraphics(frameworkHandle)
            obj.handles.framework = frameworkHandle;
            obj.params            = frameworkHandle.params;           
        end
        
        function CreateControls(obj)
           if obj.params.simulator.showSimulation

               
               
            obj.handles.graphical.mainFigure = figure('Units','normalized',...
                                                      'Tag','mainFigure',...
                                                      'ToolBar','none',...
                                                      'MenuBar','none',...                                                     
                                                      'BusyAction','queue');
                                                      
                                                  
            obj.handles.graphical.mainAxes  = axes('Units','normalized',...
                                                   'Parent',obj.handles.graphical.mainFigure,...
                                                   'NextPlot','ReplaceChildren',...
                                                   'Color',[0.2 0.1 0.23],...
                                                   'Box','on',...
                                                   'XTick',[],...
                                                   'YTick',[],...
                                                   'ZTick',[]);
             % Show domain 
               obj.params.domain.parentAxes     = obj.handles.graphical.mainAxes;
               obj.handles.classes.domain.ConstructDomain;

               
           cameratoolbar
           
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
               'Callback',@Step);
           
                obj.params.domain.parentAxes = obj.handles.graphical.mainAxes;   
              

          
            end
        end
        
        function ShowSimulation(obj)% move to Plotter class
         if obj.params.simulator.showSimulation 
            obj.ShowDomain
            obj.PlotChains 
         end
        end
                
        function PlotChains(obj)% move to Plotter class
%             if obj.params.simulator.plotChains
                for cIdx = 1:obj.params.simulator.numChains
                    p = obj.handles.classes.rouse(cIdx).position.cur;
                    set(obj.handles.graphical.chain(cIdx).beads,...
                        'XData',p(:,1),...
                        'YData',p(:,2),...
                        'ZData',p(:,3),...
                        'Visible','on');
                  
                    % plot the connectors between beads ( plot only non
                    % trivial connections)
%                     if obj.params.chain.springForce || obj.params.chain.bendingElasticityForce || obj.params.chain.lennardJonesForce
                        % Plot connectors for non trivial (non consecutives
                        % beads) connections 
%                     cm = obj.handles.classes.rouse(cIdx).connectionMap.indices.in.list;
%                     for mIdx = 1:size(cm,1);
%                         lineData.x = [obj.handles.classes.rouse(cIdx).position.cur(cm(mIdx,1)),obj.handles.classes.rouse(cIdx).position.cur(cm(mIdx,2))];
%                         lineData.y = [obj.handles.classes.rouse(cIdx).position.cur(cm(mIdx,1)),obj.handles.classes.rouse(cIdx).position.cur(cm(mIdx,2))];
%                         lineData.z = [obj.handles.classes.rouse(cIdx).position.cur(cm(mIdx,1)),obj.handles.classes.rouse(cIdx).position.cur(cm(mIdx,2))];  
%                         set(obj.handles.graphical.chain(cIdx).connectors(mIdx),...
%                             'XData',lineData.x,...
%                             'YData',lineData.y,...
%                             'ZData',lineData.z,...
%                             'Color','w',...
%                             'Tag','connector',...
%                             'Visible','on');     
%                     end                    
                    
%                     end
                    
                end
                drawnow

        end
        
        function ShowDomain(obj)% move to Plotter class
%             if obj.params.simulator.showDomain
                obj.handles.classes.domain.ShowDomain()
%                 obj.SetAxesLimits       
%             end
        end
        
        function Step(obj)
            obj.handles.framework.Step
        end
        
        function StartStop(obj,buttonHandle,varargin)% move to Plotter Class
            if obj.runSimulation
                set(buttonHandle,'String','Continue')
                obj.runSimulation = false;
            else
                set(buttonHandle,'String','Stop')
                obj.runSimulation = true;
                obj.fhandles.framework.Run
            end
        end
        
        function InitializeChainGraphics(obj)
            % set graphics
           l = linspace(0.4,0.9,obj.params.simulator.numChains);
          
           for cIdx = 1:obj.params.simulator.numChains
               r = zeros(1,3);
               for rIdx = 1:3
                k       = randperm(obj.params.simulator.numChains);
                r(rIdx) = k(1);
               end
               obj.chainColors(cIdx,1:3) = [l(r(1)),l(r(2)),l(r(3))];
           end
           
        % Create chains 
  %            cParams = obj.params.chain;
            for cIdx = 1:obj.params.simulator.numChains
                % Initialize class              
%                 obj.handles.classes.rouse(cIdx) = Rouse(cParams(cIdx));
                % make sure the chain is inside the domain
%                 obj.handles.classes.rouse(cIdx).SetInitialChainPosition(obj.handles.classes.domain);
                
                % Initialize chain graphics 
                if obj.params.simulator.showSimulation
                if ~(obj.params.chain(cIdx).springForce) && ~obj.params.chain(cIdx).bendingElasticityForce
                    lineStyle = 'none';
                else
                    lineStyle = '-';
                end
                
                [p,~] = obj.objectManager.GetPosition(cIdx);% obj.handles.classes.rouse(cIdx).position.prev;
                
                obj.handles.graphical.chain(cIdx).beads = line(...
                    'XData',p{:,1},...                                   
                    'YData',p{:,2},...
                    'ZData',p{:,3},...
                    'Marker','o',...
                    'LineStyle',lineStyle,...
                    'Color','w',...
                    'Color',obj.chainColors(cIdx,:),...
                    'Tag','chain',...
                    'MarkerSize',6,...
                    'MarkerFaceColor',obj.chainColors(cIdx,:),...
                    'Parent',obj.handles.graphical.mainAxes,...
                    'Visible','off');
                
                %TODO: fix such that the if loop is discarded
                % Plot connectors
                if any([obj.params.chain(cIdx).springForce,...
                        obj.params.chain(cIdx).bendingElasticityForce,...
                        obj.params.domain.lennardJonesForce])
                                                                        
                   cm = obj.objectManager.GetConnectionMap(cIdx);% obj.handles.classes.rouse(cIdx).connectionMap.indices.in.list;
                   cm = cm{cIdx}.indices.in.list;
                    for mIdx = 1:size(cm,1);
                        lineData.x = [obj.handles.classes.rouse(cIdx).position.cur(cm(mIdx,1)),...
                                      obj.handles.classes.rouse(cIdx).position.cur(cm(mIdx,2))];
                        lineData.y = [obj.handles.classes.rouse(cIdx).position.cur(cm(mIdx,1)),...
                                      obj.handles.classes.rouse(cIdx).position.cur(cm(mIdx,2))];
                        lineData.z = [obj.handles.classes.rouse(cIdx).position.cur(cm(mIdx,1)),...
                                      obj.handles.classes.rouse(cIdx).position.cur(cm(mIdx,2))];  
                          obj.handles.graphical.chain(cIdx).connectors(mIdx)= line(...
                            'XData',lineData.x,...
                            'YData',lineData.y,...
                            'ZData',lineData.z,...
                            'Color',obj.chainColors(cIdx,:),...
                            'Parent',obj.handles.graphical.mainAxes,...
                            'Tag','connector');     
                    end

                end
                end
            end 
        end
        
    end
end