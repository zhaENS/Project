function MoveHistonesOnChain
% This test function moves histones on a Rouse chain
close all
profile on
numSteps = 1000;
%
% % create chain and domain and register them in the ObjectManager

% Initialize simulator framework parameters
simulatorParams = SimulationFrameworkParams('showSimulation',true,'numSteps',1,'dt',0.01);

% create a spherical domain
% assign a force to the domain
sphereForces = ForceManagerParams('lennardJonesForce',false,'diffusionForce',true,'diffusionConst',0.1,...
                                  'LJPotentialWidth',0.1,'LJPotentialDepth',0.1,'dt',simulatorParams.simulator.dt);
dp(1)        = DomainHandlerParams('domainShape','sphere','forceParams',sphereForces,...
    'domainWidth',3,'dimension',simulatorParams.simulator.dimension);


% create a cylindrical Beam as a domain
cylinderForces = ForceManagerParams('diffusionForce',false,'lennardJonesForce',false,'morseForce',false);
dp(2)          = DomainHandlerParams('domainShape','cylinder','reflectionType','off','domainWidth',0.2,...
                                     'domainHeight', 25,'forceParams',cylinderForces);


% % create a chain
chainForces = ForceManagerParams('dt',simulatorParams.simulator.dt,'springForce',true,...
    'bendingElasticityForce',false,'bendingConst',1,'springConst',0.5,'minParticleEqDistance',0);
cp          = ChainParams('numBeads',100,'initializeInDomain',1,'forceParams',chainForces);
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
histoneForce = ForceManagerParams('dt',simulatorParams.simulator.dt,'diffusionConst',0.01,...
           'lennardJonesForce',false,'diffusionForce',true,'LJPotentialWidth',0.1,'LJPotentialDepth',0.1);
           
histoneParams = HistoneParams('numHistones',10,'forceParams',histoneForce);
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
    pAxes   = axes('Parent',pFigure);
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
       'MarkerFaceColor','k',...
       'markerSize',7,...
       'Parent',pAxes,...
       'LineStyle','-');
        
end
r.Run;% run initial simulator step

for sIdx = 1:numSteps
    % advance one simulation step
    r.Step;
    [~,chainPos] = r.objectManager.GetMembersPosition(1);
    chainPos     = chainPos{1};
    % move the histones
    h.Step(chainPos,simulatorParams.simulator.dt);
    
    % update histone graphics
    if simulatorParams.simulator.showSimulation
        
        set(histHandle,'XData',h.curPos(:,1),'YData',h.curPos(:,2),'ZData',h.curPos(:,3),'Parent',mAxes)
        
        % plot projected histone
        set(pHistHandle,'XData',h.curPos(:,1),'YData',h.curPos(:,2),'Parent',pAxes);
        % plot projected polymer 
        set(pPolyHandle,'XData',chainPos(:,1),'YData', chainPos(:,2),'Parent',pAxes)
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
        % update the position of the histones
        %             h.UpdateHistonePositionOnChain(chainPos)  % update the histone position on the new chain position
        % %             sIdx
    end        
end
profile viewer
end