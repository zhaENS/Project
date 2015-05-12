function MoveHistonesOnChain
% This test function moves histones on a Rouse chain 
close all 
numSteps = 1000;
% 
% % create chain and domain and register them in the ObjectManager
% % create a chain 
cp            = ChainParams('numBeads',100,'dt',0.01);
% dp            = DomainHandlerParams;
% domain        = DomainHandler(dp);
% objectManager = ObjectManager(cp);
% objectManager.InitializeObjects(domain)
% domainForceParams = dp.forceParams;
% figure, 

params = SimulationFrameworkParams(cp);

% SimulatorParams
% profile on 
r = RouseSimulatorFramework(params);
[~,initialChainPosition] = r.objectManager.GetMembersPosition(1);
initialChainPosition     = initialChainPosition{1};
% Initialize histones
histoneParams.dt             = r.params.simulator.dt;
histoneParams.numHistones    = 10;
histoneParams.diffusionConst = 0.1;
h                            = Histone(histoneParams, initialChainPosition);

% Test with one histone 
mAxes = r.simulationGraphics.handles.graphical.mainAxes;
% add initial histone position 
histHandle = line('XData',h.curPos(:,1),'YData',h.curPos(:,2),'ZData',h.curPos(:,3),'marker','o','MarkerFaceColor','y','MarkerSize',10,'Parent',mAxes,'LineStyle','none');

    r.Run;
    for sIdx = 1:numSteps
        r.Step;
        [~,chainPos] = r.objectManager.GetMembersPosition(1);
        chainPos = chainPos{1};
        h.Step(chainPos);
        
        % add the histone to the graphics         
        set(histHandle,'XData',h.curPos(:,1),'YData',h.curPos(:,2),'ZData',h.curPos(:,3))        
    end    
end