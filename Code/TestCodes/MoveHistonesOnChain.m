function MoveHistonesOnChain
% This test function moves histones on a Rouse chain
close all
% profile on
numSteps = 1000;
%
% % create chain and domain and register them in the ObjectManager

% Initialize simulator framework parameters
simulatorParams = SimulationFrameworkParams('showSimulation',true,'numSteps',1,'dt',0.1);

% create a spherical domain
% assign a force to the domain
sphereForces = ForceManagerParams('lennardJonesForce',false,'diffusionForce',false,'diffusionConst',0.001,...
                                  'mechanicalForce',true,'mechanicalForceDirection','out',...
                                  'mechanicalForceCenter',[0 0 0],'mechanicalForceMagnitude',0.01,...
                                  'LJPotentialWidth',0.1,'LJPotentialDepth',0.1,'dt',simulatorParams.simulator.dt);
dp(1)        = DomainHandlerParams('domainShape','sphere','forceParams',sphereForces,...
                                   'domainWidth',0.1,'dimension',simulatorParams.simulator.dimension);


% create a cylindrical Beam as a domain
cylinderForces = ForceManagerParams('diffusionForce',false,'lennardJonesForce',false,'morseForce',false);
dp(2)          = DomainHandlerParams('domainShape','cylinder','reflectionType','off','domainWidth',0.01,...
                                     'domainHeight', 10,'forceParams',cylinderForces);

% % create a chain
chainForces = ForceManagerParams('dt',simulatorParams.simulator.dt,'springForce',true,...
                    'bendingElasticityForce',false,'bendingConst',1,'springConst',1,'minParticleEqDistance',0);
cp          = ChainParams('numBeads',500,'initializeInDomain',1,'forceParams',chainForces,'b',0.01);
% cp(2)     = ChainParams('numBeads',100,'initializeInDomain',1,'forceParams',chainForces);

% register the object parameters in the simulator framework
simulatorParams.SetDomainParams(dp);
simulatorParams.SetChainParams(cp);

% initialize simulator framework
r = RouseSimulatorFramework(simulatorParams);

% get the chain position to initialize the histone on it
[~,initialChainPosition] = r.objectManager.GetMembersPosition(1);
initialChainPosition     = initialChainPosition{1};

% Initialize histones with the chain position
histoneForce = ForceManagerParams('dt',simulatorParams.simulator.dt,'diffusionConst',0.01,'mechanicalForce',false,...
                       'mechanicalForceDirection','out','mechanicalForceMagnitude',0.3,'mechanicalForceCenter',[0 0 0],...
                       'lennardJonesForce',false,'diffusionForce',false,'LJPotentialWidth',0.01,'LJPotentialDepth',0.01);
           
histoneParams = HistoneParams('numHistones',200,'forceParams',histoneForce);
h             = Histone(histoneParams);

h.Initialize(initialChainPosition);

% get axes
if simulatorParams.simulator.showSimulation
    mAxes = r.simulationGraphics.handles.graphical.mainAxes;
   
    % initialize histone graphics %TODO: incorporate histone graphics in the simulationGraphics class
    histHandle = line('XData',h.curPos(:,1),...
                      'YData',h.curPos(:,2),...
                     'ZData',h.curPos(:,3),...
                     'marker','o',...
                     'MarkerFaceColor','y',...
                     'MarkerSize',10,...
                     'Parent',mAxes,...
                     'LineStyle','none');
    daspect(mAxes,[1 1 1])
    % create figure for the projection in the x-y plane 
    pFigure = figure;
    pAxes   = axes('Parent',pFigure,'XLim',get(mAxes,'XLim'), 'YLim',get(mAxes,'YLim'),'FontSize',30);
    % projected histone 
    pHistHandle = line('XData',h.curPos(:,1),...
                     'YData',h.curPos(:,2),...
                     'Marker','o',...
                     'MarkerFaceColor','y',...
                     'MarkerSize',10,...
                     'Parent',pAxes,...
                     'LineStyle','none');
   % projected Polymer 
   pPolyHandle = line('XData',initialChainPosition(:,1),...
                      'YDAta',initialChainPosition(:,2),...
                      'Marker','o',...
                      'MarkerFaceColor','none',...
                      'MarkerEdgeColor','b',...
                      'markerSize',7,...
                      'Parent',pAxes,...
                      'LineStyle','-');
  
   % create density axes 
   dFig = figure;
   dAxes = axes('Parent',dFig,'NextPlot','add');
   xlabel(dAxes,'Time','FontSize',40);
   ylabel(dAxes,'Density','FontSize',40)
    % insert roi rectangle 
   rectX      = -0.1; % buttom left corner
   rectY      = -0.1; % buttom left corner
   rectWidth  = 0.2;
   rectHeight = 0.2;
   % insert the projection to the 3d axes
    patch([rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)],...
        'r', 'Parent',mAxes, 'FaceAlpha',0.5);
   % insert the projection to the projection axes
    patch([rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)],...
        'r', 'Parent',pAxes, 'FaceAlpha',0.5);     
end
r.Run;% run initial simulator step

for sIdx = 1:numSteps
    % Advance one simulation step
    r.Step;
    [~,chainPos] = r.objectManager.GetMembersPosition(1);
    chainPos     = chainPos{1};
    % move the histones
    h.Step(chainPos,simulatorParams.simulator.dt);% update current position 
       
    if simulatorParams.simulator.showSimulation
        % update histone graphics
        set(histHandle,'XData',h.curPos(:,1),'YData',h.curPos(:,2),'ZData',h.curPos(:,3));        
        % plot projected histone
        set(pHistHandle,'XData',h.curPos(:,1),'YData',h.curPos(:,2));
        % plot projected polymer 
        set(pPolyHandle,'XData',chainPos(:,1),'YData', chainPos(:,2))
        drawnow
    end
    
    % Apply bending elasticity forces for beads inside the beam
    inBeam = r.handles.classes.domain.InDomain(chainPos,2);
    if any(inBeam)
        bendingElasticityConst = 1/simulatorParams.simulator.dt;
        connectivityMat        = r.objectManager.GetConnectivityMapAsOne(1);
        bForce                 = ForceManager.GetBendingElasticityForce(true,chainPos,connectivityMat,bendingElasticityConst ,[]);
        
        % zero out forces outside the beam
        bForce(~inBeam,:) = 0;
        
        % update the position of the chain
        chainPos(inBeam,:) = chainPos(inBeam,:) + bForce(inBeam,:)*simulatorParams.simulator.dt;

        % assign the new position to the chain
        r.objectManager.DealCurrentPosition(1,chainPos);

    end   
    % calculate the density of DNA and histones in the roi
    % find all beads in the roi     
%     chainLength = cumsum(sqrt(sum(chainPos.^2,2)));
%     chainLength = chainLength(end);
%     beadsIn     = (chainPos(:,1)<=1 & chainPos(:,1)>=-1 & chainPos(:,2)>=-1 & chainPos(:,2)<=1);
%     beadsOut    = ~beadsIn;
    
    
    % calculate the histoneDensity        
    histoneDensity =  sum((h.curPos(:,1)<=(rectX+rectWidth) & h.curPos(:,1)>=rectX &...
                           h.curPos(:,2)<=(rectY+rectHeight) & h.curPos(:,2)>=rectY))/h.params.numHistones;
    line('XData',r.simulationData.step,'YData',histoneDensity,'Marker','.','Parent',dAxes)
end
% profile viewer
end