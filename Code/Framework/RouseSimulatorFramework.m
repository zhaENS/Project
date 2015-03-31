classdef RouseSimulatorFramework<handle
    %TODO: debug 
    %TODO: add plotter class
    %TODO: add presets for simulations,(load/save)
    %TODO: move the initialization of the chains to allow initialization on
    %the fly
    %TODO: add a listener that draws emails from the server, to allow
    %controlling the simulations from a distance
    %TODO: save noise terms seeds to allow resimulations
    %TODO: save parameters and a short explanation with the results
    %TODO: parameters should have their own class for each class
    %TODO: remove the option to plot online. plotting is always offline.
    %remove parentAxes property from domainHandlerParams
    %TODO: chain parameters should be initialized before the chains according to the number specified by the simulator params 
    %TODO: add class for each force. input: position, connectivityMat,
    %distances. output: bead displacement
    %TODO: 
    properties
        handles        
        beadDistance
        runSimulation
        params
        simulationData  = struct('step',[],'stepTime',[]);      
        batchRound      = 0;
        simulationRound = 0; 
    end
    
    properties (Access=private)
        chainColors 
        chainsCenterOfMass
        recipe = struct('PreSimulationBatchActions','',...
                        'PreRunActions','',...
                        'PostRunActions','',...
                        'PostSimulationBatchActions','');
    end
    
    events
        run        
    end
    
    methods 
        function obj = RouseSimulatorFramework(frameworkParams)
                                
           % Set initial parameters
           obj.OrganizeParams(frameworkParams)
           % Create controls
           obj.CreateControls 
           % Initizlie chains and domain 
           obj.InitializeClasses
           % set recipe files 
           obj.ReadRecipeFile           
        end
        
        function OrganizeParams(obj,frameworkParams)
            obj.params                             = frameworkParams;            
            obj.params.chain.dimension             = obj.params.simulator.dimension;
            obj.params.domain.showDomain           = obj.params.simulator.showSimulation;  
            obj.params.domain.dimension            = obj.params.simulator.dimension;
            obj.params.dataRecorder.recipeFileName = obj.params.simulator.recipeFileName;
            obj.params.dataRecorder.encounterDist  = obj.params.simulator.encounterDist;
        end
                
        function SetInputParams(obj,varargin)
            if mod(numel(varargin),2)~=0
                error('value pair input must come in pairs')
            end
            for pIdx = 1:2:numel(varargin);
                obj.params.(varargin{pIdx}) = varargin{pIdx+1};
            end
        end
                
        function SetChainParams(obj)%obsolete
                       % set chain parameters 
           for cIdx = 1:obj.params.numChains
               obj.params.chain{cIdx}    = ...
                         {'dt',obj.params.dt,... 
                          'numBeads',obj.params.numBeads,...
                          'dimension',obj.params.dimension,...
                          'diffusionConst',obj.params.diffusionConst,...   
                          'minBeadDist',obj.params.minBeadDist,...
                          'numSteps',obj.params.numSteps,...
                          'showSimulation',false,...
                          'debug',false,...
                          'showBeadsDist',false,...
                          'showEnd2EndDist',false,...
                          'fixedBeadNum',obj.params.fixedBeadNum,...
                          'b',obj.params.b,...
                          'connectMonomers',obj.params.connectMonomers,...
                          'bendingElasticityForce',obj.params.bendingElasticityForce,...  
                          'bendingConst',obj.params.bendingConst,...
                          'lennardJonesForce',obj.params.lennardJonesForce,...
                          'springForce',obj.params.springForce,...
                          'LJPotentialDepth',obj.params.LJPotentialDepth,...
                          'LJPotentialWidth',obj.params.LJPotentialWidth};                                                                       
           end 
        end
        
        function InitializeClasses(obj)
            
            obj.CreateDomain
            
            obj.InitializeChains                                          
%             obj.InitializeDistributionHandler
            obj.InitializeDataRecorder;
        end
        
        function InitializeChains(obj)
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
            for cIdx = 1:obj.params.simulator.numChains
                % Initialize class
                cParams                         = obj.params.chain;% should be expanded to include parameters for each chain
                obj.handles.classes.rouse(cIdx) = Rouse(cParams);
                % make sure the chain is inside the domain
                obj.handles.classes.rouse(cIdx).SetInitialChainPosition(obj.handles.classes.domain);
                if ~(obj.params.chain.springForce) && ~obj.params.chain.bendingElasticityForce
                    lineStyle = 'none';
                else
                    lineStyle = '-';
                end
                % Initialize chain graphics
                p = obj.handles.classes.rouse(cIdx).positions.beads.prev;
                obj.handles.graphical.chain(cIdx).beads = line(...
                    'XData',p.x,...                                   
                    'YData',p.y,...
                    'ZData',p.z,...
                    'Marker','o',...
                    'LineStyle',lineStyle,...
                    'Color','w',...
                    'Color',obj.chainColors(cIdx,:),...
                    'Tag','chain',...
                    'MarkerSize',6,...
                    'MarkerFaceColor',obj.chainColors(cIdx,:),...
                    'Parent',obj.handles.graphical.mainAxes,...
                    'Visible','off');
                
                
                % Plot connectors
                if obj.params.chain.springForce || obj.params.chain.bendingElasticityForce
                                                                        
                   cm = obj.handles.classes.rouse(cIdx).connectionMap.indices.in.list;
                    for mIdx = 1:size(cm,1);
                        lineData.x = [obj.handles.classes.rouse(cIdx).positions.beads.cur.x(cm(mIdx,1)),...
                                      obj.handles.classes.rouse(cIdx).positions.beads.cur.x(cm(mIdx,2))];
                        lineData.y = [obj.handles.classes.rouse(cIdx).positions.beads.cur.y(cm(mIdx,1)),...
                                      obj.handles.classes.rouse(cIdx).positions.beads.cur.y(cm(mIdx,2))];
                        lineData.z = [obj.handles.classes.rouse(cIdx).positions.beads.cur.z(cm(mIdx,1)),...
                                      obj.handles.classes.rouse(cIdx).positions.beads.cur.z(cm(mIdx,2))];  
                          obj.handles.graphical.chain(cIdx).connectors(mIdx)= line(...
                            'XData',lineData.x,...
                            'YData',lineData.y,...
                            'ZData',lineData.z,...
                            'Color',obj.chainColors(cIdx,:),...
                            'Parent',obj.handles.graphical.mainAxes,...
                            'Tag','connector');     
                    end
%                 end
                end
            end                                 
           % get chains center of mass- for the plotting 
           obj.GetChainsCenterOfMass 
        end
        
        function ReadRecipeFile(obj)
%             try
         t = fileread([obj.params.simulator.recipeFileName,'.rcp']);
%             catch
                % if the file does not exist, read the default file 
%             end
         % search for the function marker 
         [funcStartPos1,funcStartPos2] = regexp(t,'<func>');
         [funcEndPos1,funcEndPos2]     = regexp(t,'</func>');
         for fIdx = 1:numel(funcStartPos1)-1
             % sort the functions into categories
             functionName = strtrim(t(funcStartPos2(fIdx)+1:funcEndPos1(fIdx)-1));
             obj.recipe.(functionName) = t(funcEndPos2(fIdx)+1:funcStartPos1(fIdx+1)-1);
         end
             functionName = strtrim(t(funcStartPos2(end)+1:funcEndPos1(end)-1));
             obj.recipe.(functionName) = t(funcEndPos2(end)+1:end);     
        end
               
        function InitializeDataRecorder(obj)
            obj.handles.classes.dataRecorder = SimulationDataRecorder(obj.params.dataRecorder);
        end
        
        function PreSimulationBatchActions(obj)
            % Actios performed before a simulation batch      
            obj.batchRound      = obj.batchRound+1;
            obj.simulationRound = 0;
            obj.handles.classes.dataRecorder.CreateNewSimulationBatch
            
            % evaluate the recipe file at the PreSimulationBatchActions
            eval(obj.recipe.PreSimulationBatchActions)
        end
        
        function PreRunActions(obj)
            % Actions called before running the simulation        
            obj.simulationRound = obj.simulationRound+1;
            obj.runSimulation   = true;            
            
            eval(obj.recipe.PreRunActions);
            
            obj.handles.classes.dataRecorder.NewSimulation(obj.handles.classes.rouse,obj.params);                        
            obj.handles.classes.dataRecorder.SetSimulationStartTime;
            obj.simulationData(obj.batchRound,obj.simulationRound).step     = 1;
            obj.simulationData(obj.batchRound,obj.simulationRound).stepTime = 0;
        end
        
        function Run(obj)
            % Run the simulation according to the specified simulationType
            % parameter
           
            for bIdx = 1:obj.params.simulator.numSimulationBatches
                % perform actions before each simulation batch
                obj.PreSimulationBatchActions
                for sIdx = 1:obj.params.simulator.numSimulations
                    % perform actions before each simulation
                    obj.PreRunActions                    
                    while obj.runSimulation
                        % perform actions before the current step
                        obj.PreStepActions                                               
                        obj.NextStep;% advance one step of the polymer chain                           
                        % perform action post the current step 
                        obj.PostStepActions                                              
                    end
                    % perform actions post simulation
                    obj.PostRunActions
                end
                % perform actions post simulation batch
                obj.PostSimulationBatchActions
            end
           
        end               
        
        function PreStepActions(obj)
            eval(obj.recipe.PreStepActions);
        end
         
        function PostStepActions(obj)
            eval(obj.recipe.PostStepActions);
              obj.Record % record data related to the Rouse polymer                        
              % raise the Stop flag if the number of steps is more than the allowed
              obj.runSimulation = obj.runSimulation && ...
             (obj.simulationData(obj.batchRound,obj.simulationRound).step<obj.params.simulator.numSteps);
        end
        
        function PostRunActions(obj)
            % Calculate distributions related to the Rouse polymer
            % chain
            
            eval(obj.recipe.PostRunActions);
            
            obj.handles.classes.dataRecorder.SetSimulationEndTime;
            obj.handles.classes.dataRecorder.SaveResults;% save results
            
            if obj.params.simulator.notifyByEmail
                if mod(obj.simulationRound,obj.params.simulator.notifyCycleLength)==0
                    msg = sprintf('%s%s%s\n%s%s%s\n%s%s\n%s','simulation batch',num2str(obj.batchRound),', ',...
                        'simulation ', num2str(obj.simulationRound),' is done. ',...
                        'Simulation Time: ',num2str(obj.handles.classes.dataRecorder.simulationData(obj.simulationRound).simulationTime),...
                        'Simulating the distance between monomenrs at relexation time, each time bead 1 is connected to bead j, where j=3..64');
                    try
                        SendMail('ofirshukron','Simulation Status',msg);
                    catch
                    end
                end
            end
        end
        
        function PostSimulationBatchActions(obj)
            % actions performed post simulation batch
            eval(obj.recipe.PostSimulationBatchActions);            
            obj.handles.classes.dataRecorder.ClearAllSimulationData;            
%             obj.simulationRound = 0;               
        end
        
        function CreateControls(obj)% move to Plotter class
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
               'Callback',@obj.NextStep);
                      
            end
        end
        
        function CreateDomain(obj)
            
             % parameters for the domain are set in the organizeParams
             obj.params.domain.parentAxes = obj.handles.graphical.mainAxes;
             obj.handles.classes.domain  = DomainHandler(obj.params.domain);
             
            if obj.params.domain.showDomain
                % move these commands to future Plotter class               
               obj.params.domain.parentAxes     = obj.handles.graphical.mainAxes;
               obj.handles.classes.domain.ConstructDomain;
               set(obj.handles.graphical.mainAxes,'NextPlot','Add')                 
               set(obj.handles.graphical.mainAxes,'NextPlot','ReplaceChildren')
            end

        end
        
        function StartStop(obj,buttonHandle,varargin)% move to Plotter Class
            if obj.runSimulation 
                set(buttonHandle,'String','Continue')
                obj.runSimulation = false;
            else 
                set(buttonHandle,'String','Stop')
                obj.runSimulation = true;
                obj.Run
            end
        end
        
        function NextStep(obj,varargin)
            % Next simulation step 
           
            for cIdx = 1:obj.params.simulator.numChains 
                obj.handles.classes.rouse(cIdx).Next % advance one step for all chains 
                obj.Reflect(cIdx) % reflection with the domain             
%                 obj.handles.classes.rouse(cIdx).GetForces;
                obj.handles.classes.rouse(cIdx).SetPrevBeadPosition;

            end            
            obj.ShowSimulation            
            obj.simulationData(obj.batchRound,obj.simulationRound).step = ...
                obj.simulationData(obj.batchRound,obj.simulationRound).step+1;

        end
        
        function Reflect(obj,chainIdx)
            % Reflect the beads
            dimName  = {'x','y','z'};            
            curPos   = obj.handles.classes.rouse(chainIdx).positions.beads.cur; 
            prevPos  = obj.handles.classes.rouse(chainIdx).positions.beads.prev; 
            
            inIdx    = obj.handles.classes.domain.InDomain(curPos);
            outIdx   = find(~inIdx);
            for oIdx = 1:numel(outIdx)
                for dIdx = 1:3
                    beadCur.(dimName{dIdx})  = curPos.(dimName{dIdx})(outIdx(oIdx));
                    beadPrev.(dimName{dIdx}) = prevPos.(dimName{dIdx})(outIdx(oIdx));
                end
                [~,np] = obj.handles.classes.domain.Reflect(beadPrev,beadCur);
                for dIdx = 1:obj.params.simulator.dimension
                    obj.handles.classes.rouse(chainIdx).positions.beads.cur.(dimName{dIdx})(outIdx(oIdx))  = np.(dimName{dIdx});                    
                end                 
            end  
        end        
        
        function Record(obj)
            % Add samples to the recorded distributions in
            % distributionHandler
            if obj.params.simulator.recordData
%                obj.handles.classes.distributionHandler.Add% add to the frequency array of the distributionHandler
               obj.handles.classes.dataRecorder.Add;
            end
        end
        
        function HighlightCloseBeads(obj,chainInd1,chainInd2)% move to Plotter class
            for bIdx = 1:numel(chainInd1)
                beads.x(2*bIdx-1) = obj.handles.classes.rouse(chainInd1(bIdx)).positions.beads.cur.x(end,end);
                beads.x(2*bIdx)   = obj.handles.classes.rouse(chainInd2(bIdx)).positions.beads.cur.x(end,end);
                beads.y(2*bIdx-1) = obj.handles.classes.rouse(chainInd1(bIdx)).positions.beads.cur.y(end,end);
                beads.y(2*bIdx)   = obj.handles.classes.rouse(chainInd2(bIdx)).positions.beads.cur.y(end,end);
                beads.z(2*bIdx-1) = obj.handles.classes.rouse(chainInd1(bIdx)).positions.beads.cur.z(end,end);
                beads.z(2*bIdx)   = obj.handles.classes.rouse(chainInd2(bIdx)).positions.beads.cur.z(end,end);
            end  
            set(obj.handles.graphical.mainAxes,'NextPlot','Add')
                line('XData',beads.x,...
                    'YData',beads.y,...
                    'ZData',beads.z,...
                    'Parent',obj.handles.graphical.mainAxes,...
                    'MarkerEdgeColor','r',...
                    'Marker','o',...
                    'LineStyle','-',...
                    'Color','r',...
                    'Tag','closeBeads',...
                    'MarkerSize',9);                    
                
            set(obj.handles.graphical.mainAxes,'NextPlot','ReplaceChildren')

            
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
                    p = obj.handles.classes.rouse(cIdx).positions.beads.cur;
                    set(obj.handles.graphical.chain(cIdx).beads,...
                        'XData',p.x,...
                        'YData',p.y,...
                        'ZData',p.z,...
                        'Visible','on');
                    
                    % plot the connectors between beads ( plot only non
                    % trivial connections)
                    if obj.params.chain.springForce || obj.params.chain.bendingElasticityForce
                        % Plot connectors for non trivial (non consecutives
                        % beads) connections 
                    cm = obj.handles.classes.rouse(cIdx).connectionMap.indices.in.list;
                    for mIdx = 1:size(cm,1);
                        lineData.x = [obj.handles.classes.rouse(cIdx).positions.beads.cur.x(cm(mIdx,1)),obj.handles.classes.rouse(cIdx).positions.beads.cur.x(cm(mIdx,2))];
                        lineData.y = [obj.handles.classes.rouse(cIdx).positions.beads.cur.y(cm(mIdx,1)),obj.handles.classes.rouse(cIdx).positions.beads.cur.y(cm(mIdx,2))];
                        lineData.z = [obj.handles.classes.rouse(cIdx).positions.beads.cur.z(cm(mIdx,1)),obj.handles.classes.rouse(cIdx).positions.beads.cur.z(cm(mIdx,2))];  
                        set(obj.handles.graphical.chain(cIdx).connectors(mIdx),...
                            'XData',lineData.x,...
                            'YData',lineData.y,...
                            'ZData',lineData.z,...
                            'Color','w',...
                            'Tag','connector',...
                            'Visible','on');     
                    end                    
                    
                    end
                    
                end
                drawnow
%                 start(obj.handles.classes.timer);
%                 stop(obj.handles.classes.timer);
                
%             end
        end
        
        function ShowDomain(obj)% move to Plotter class
%             if obj.params.simulator.showDomain
                obj.handles.classes.domain.ShowDomain()
%                 obj.SetAxesLimits       
%             end
        end
        
        function SetAxesLimits(obj)% move to Plotter class
            if strcmpi(obj.params.domainShape,'none')
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
        
        function GetChainsCenterOfMass(obj)% move to chain class
           % get the center of mass of ALL chains  
%            ccm(1:obj.params.simulator.numChains) = struct('x',[],'y',[],'z',[]);
%            
%             for cIdx = 1:obj.params.simulator.numChains              
%               ccm(cIdx) = obj.handles.classes.rouse(cIdx).positions.beads.cm;             
%             end
%            
%            cm.x = sum([ccm.x])/obj.params.simulator.numChains;
%            cm.y = sum([ccm.y])/obj.params.simulator.numChains;
%            cm.z = sum([ccm.z])/obj.params.simulator.numChains;
           obj.chainsCenterOfMass = zeros(1,obj.params.chain.dimension);%cm;
        end
    end
       
end
    