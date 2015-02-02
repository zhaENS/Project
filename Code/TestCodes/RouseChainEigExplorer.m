classdef RouseChainEigExplorer<handle
    % this class creates a Gui to allow the visualization of the effect of
    % changing the Rouse eigenvalues on the connectivity (Kirchhoff matrix)
    % of the chain. 
    % A linear chain is first formed. R=V'EV. The eigenvectors V and eigen vealues E are computed and the
    % user is allowed to change the E values to form a new matrix E_c . It later displayes the
    % resuling Kirchhoff matrix by R= VE_cV';
    % note that the eigenvetors do not change. 
    
properties
    handles
    matrix% the rouse matrix
    eigenvalues
    eigenvectors
end

methods
    function obj = RouseChainEigExplorer
        obj.matrix = RouseMatrix(20);
        [obj.eigenvectors, obj.eigenvalues] = eig(obj.matrix);
        
         obj.handles.mainFigure = figure('Units','norm','MenuBar','none','ToolBar','none');
         % display panel
         obj.handles.displayPanel= uipanel('Parent',obj.handles.mainFigure,...
             'Units','norm',...             
             'Position',[0 0 0.9 1]);
         obj.handles.mainAxes   = axes('Parent',obj.handles.displayPanel,...
             'Units','norm',...
             'Position',[0.05 0.05 0.90 0.90],...
             'NextPlot','ReplaceChildren');
           
         % controlPanel 
         obj.handles.controlPanel = uipanel('Parent',obj.handles.mainFigure,...
             'Units','norm',...
             'Position',[0.9 0 0.1 1]);
             
         % slider 
         obj.handles.slider  = uicontrol('Parent',obj.handles.controlPanel,...
             'Units','norm',...
             'Position',[0 0.3 0.7, 0.5],...
             'SliderStep',[0.001,0.01],...
             'Callback',@obj.SliderCallback,...
             'Style','slider',...
             'Min', 0,...
             'Max',2*(size(obj.matrix,1)-1));
         
         % eig index
         eigString = num2str(1);
         for bIdx = 2:size(obj.matrix,1)
           eigString = [eigString '|', num2str(bIdx)];
         end
         
         obj.handles.dropDown = uicontrol('Parent',obj.handles.controlPanel,...
             'Style','popupmenu',...
             'Units','norm',...
             'Position',[0 0.2, 0.7,0.1],...
             'String',eigString,...             
             'Callback',@obj.ChooseEig);
         
         obj.handles.editEigVal =  uicontrol('Parent',obj.handles.controlPanel,...
             'Units','norm',...
             'Position',[0,0.1 0.7 0.1],...
             'Style','edit',...
             'String','',...
             'Callback',@obj.EditEigValue);
         
         % eigenvalue value
         obj.handles.eigVal = uicontrol('Parent',obj.handles.controlPanel,...
             'Units','norm',...
             'Position',[0,0.8 0.7 0.1],...
             'Style','text');
         
         % num beads
         obj.handles.numBeads = uicontrol('Parent',obj.handles.controlPanel,...
             'Units','norm',...
             'Position',[0,0.9 0.7 0.1],...
             'Style','edit',...
             'String',size(obj.matrix,1),...
             'Callback',@obj.NumBeadsCallback);
         
         obj.Display
    end
    
    function Display(obj)    
        c = get(obj.handles.mainAxes,'CameraPosition');
        surf(obj.handles.mainAxes,obj.matrix);%,[],'Parent',obj.handles.mainAxes) 
        set(obj.handles.mainAxes,'CameraPosition',c)
        cameratoolbar
    end
    
    function EditEigValue(obj,varargin)
         eigInd = (get(obj.handles.dropDown,'Value'));
         newVal = str2double(get(obj.handles.editEigVal,'String'));
         set(obj.handles.slider,'Value',newVal);
         obj.eigenvalues(eigInd,eigInd)= newVal;          
         obj.matrix = obj.eigenvectors*obj.eigenvalues*obj.eigenvectors';
         obj.Display                          
    end
    
    function NumBeadsCallback(obj,varargin)
        % change the number of beads 
        numBeads = str2double(get(obj.handles.numBeads,'String'));
        obj.matrix = RouseMatrix(numBeads);
        [obj.eigenvectors, obj.eigenvalues] = eig(obj.matrix);
        obj.eigenvalues = max(obj.eigenvalues,0);
        eigString = num2str(1);
         for bIdx = 2:size(obj.matrix,1)
           eigString = [eigString '|', num2str(bIdx)];
         end
         set(obj.handles.dropDown,'String',eigString);
         set(obj.handles.dropDown,'Value',1);
         set(obj.handles.slider,'Value',obj.eigenvalues(1,1));
%          set(obj.handles.slider,'Max',numBeads-1);         
        obj.Display
    end
    
    function SliderCallback(obj,varargin)
        sliderVal = get(obj.handles.slider,'Value');
        eigInd = (get(obj.handles.dropDown,'Value'));
        obj.eigenvalues(eigInd,eigInd)= sliderVal;
        set(obj.handles.editEigVal,'String',num2str(obj.eigenvalues(eigInd,eigInd)));
        obj.matrix = obj.eigenvectors*obj.eigenvalues*obj.eigenvectors';
        obj.Display
    end
    
    function ChooseEig(obj,varargin)
        eigInd = (get(obj.handles.dropDown,'Value'));
        set(obj.handles.editEigVal,'String',num2str(obj.eigenvalues(eigInd,eigInd)));
        set(obj.handles.slider,'Value',obj.eigenvalues(eigInd,eigInd));
%         set(obj.handles.eigVal,'String',obj.eigenvalues(eigInd,eigInd));
    end
end


end
