function MoveHistonesOnChain
% This test function moves histones on a Rouse chain 
close all 
numSteps = 1000;
% 
% % create chain and domain and register them in the ObjectManager
% % create a chain 
cp     = ChainParams('numBeads',100,'dt',0.01);

params = SimulationFrameworkParams(cp);
params.simulator.numSteps = 1;% one step for initialization 

% initialize simulator framework
r = RouseSimulatorFramework(params);

[~,initialChainPosition] = r.objectManager.GetMembersPosition(1);
initialChainPosition     = initialChainPosition{1};

% Initialize histones
h                        = Histone('dt',r.params.simulator.dt,'numHistones',1,'diffusionConst',1,'ljForce',false,...
                                   'diffusionForce',true,'ljPotentialWidth',0.1,'ljPotentialDepth',0.1);
h.Initialize(initialChainPosition);

% get axes
mAxes = r.simulationGraphics.handles.graphical.mainAxes;

% initialize histone graphics
histHandle = line('XData',h.curPos(:,1),...
                  'YData',h.curPos(:,2),...
                  'ZData',h.curPos(:,3),...
                  'marker','o',...
                  'MarkerFaceColor','y',...
                  'MarkerSize',10,...
                  'Parent',mAxes,...
                  'LineStyle','none');

    r.Run;% run initial simulator step 
    
    for sIdx = 1:numSteps
        r.Step; % advance one simulation step 
        [~,chainPos] = r.objectManager.GetMembersPosition(1);
        chainPos     = chainPos{1};
        h.Step(chainPos); % move the histones
        
        % update histone graphics         
        set(histHandle,'XData',h.curPos(:,1),'YData',h.curPos(:,2),'ZData',h.curPos(:,3))        
    end    
end