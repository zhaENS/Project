function MoveHistonesOnChain
% This test function moves histones on a Rouse chain 
close all 
numSteps = 200;
% 
% % create chain and domain and register them in the ObjectManager

% create a spherical domain 
% assign a force to the domain 
sphereForces = ForceManagerParams('lennardJonesForce',false,'diffusionForce',true,'diffusionConst',1,...
                                  'LJPotentialWidth',0.1,'LJPotentialDepth',0.1);
dp(1)        = DomainHandlerParams('domainShape','sphere','forceParams',sphereForces);
                               

% create a cylindrical Beam as a domain
cylinderForces = ForceManagerParams('diffusionForce',false,'lennardJonesForce',false,'morseForce',false);                                
dp(2)          = DomainHandlerParams('domainShape','cylinder','reflectionType','off','domainWidth',1,...
                                     'domainHeight', 50,'forceParams',cylinderForces);
                                

% % create a chain 

chainForces = ForceManagerParams('dt',0.01,'springForce',true,'bendingElasticityForce',false,'bendingConst',1,'springConst',1);
cp          = ChainParams('numBeads',100,'initializeInDomain',1,'forceParams',chainForces);
% cp(2)     = ChainParams('numBeads',100,'dt',0.01,'initializeInDomain',1);

% initialize simulator framework
simulatorParams = SimulationFrameworkParams();
simulatorParams.simulator.numSteps = 1;% one step for initialization 
simulatorParams.SetDomainParams(dp);
simulatorParams.SetChainParams(cp);
r = RouseSimulatorFramework(simulatorParams);

% get the chain position to initialize the histone on it 
[~,initialChainPosition] = r.objectManager.GetMembersPosition(1);
initialChainPosition     = initialChainPosition{1};

% Initialize histones with the chain position 
histoneForce = ForceManagerParams('dt',r.params.simulator.dt,'diffusionConst',4,...
                                  'lennardJonesForce',false,'diffusionForce',true,'LJPotentialWidth',0.1,'LJPotentialDepth',0.1);    
h                        = Histone('numHistones',10,'forceParams',histoneForce);
                                  
h.Initialize(initialChainPosition);

% get axes
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
daspect([1 1 1])
    r.Run;% run initial simulator step 
    
    for sIdx = 1:numSteps
        r.Step; % advance one simulation step 
        [~,chainPos] = r.objectManager.GetMembersPosition(1);
        chainPos     = chainPos{1};
        h.Step(chainPos,simulatorParams.simulator.dt); % move the histones                
        
        % Apply bending elasticity forces for beads inside the beam
        inBeam = r.handles.classes.domain.InDomain(chainPos,2);
        if any(inBeam)            
            bendingElasticityConst = 1/simulatorParams.simulator.dt;
            bForce = ForceManager.GetBendingElasticityForce(true,chainPos,r.objectManager.connectivity,bendingElasticityConst ,[]);
            % zero out forces outside the beam 
            bForce(~inBeam,:) = 0;
            % update the position of the chain 
            chainPos(inBeam,:) = chainPos(inBeam,:) + bForce(inBeam,:)*simulatorParams.simulator.dt;
            % assign the new position to the chain 
            r.objectManager.DealCurrentPosition(1,chainPos);
            % update the position of the histones
            h.UpdateHistonePositionOnChain(chainPos)  % update the histone position on the new chain position 
            sIdx
        end
        
        % update histone graphics         
        set(histHandle,'XData',h.curPos(:,1),'YData',h.curPos(:,2),'ZData',h.curPos(:,3))        
    end    
end