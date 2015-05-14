function MoveHistonesOnChain
% This test function moves histones on a Rouse chain 
close all 
numSteps = 500;
% 
% % create chain and domain and register them in the ObjectManager
% create a domain 
dp(1)     = DomainHandlerParams('domainShape','sphere','diffusionForce',true,'diffusionConst',0.05,...
                                'lennardJonesForce',true);

% create a Beam as a domain
dp(2)     = DomainHandlerParams('domainShape','cylinder','reflectionType','off','domainWidth',3,...
                                'domainHeight', 50,'diffusionForce',false,'lennardJonesForce',false,...
                                'morseForce',false);
                                

% % create a chain 
cp     = ChainParams('numBeads',32,'dt',0.01,'initializeInDomain',1,'springForce',true,'bendingElasticityForce',false);
% cp(2)     = ChainParams('numBeads',100,'dt',0.01,'initializeInDomain',1);

params = SimulationFrameworkParams(cp,dp);
params.simulator.numSteps = 1;% one step for initialization 

% initialize simulator framework
r = RouseSimulatorFramework(params);

% get the chain position to initialize the histone on it 
[~,initialChainPosition] = r.objectManager.GetMembersPosition(1);
initialChainPosition     = initialChainPosition{1};

% Initialize histones with the chain position 
h                        = Histone('dt',r.params.simulator.dt,'numHistones',10,'diffusionConst',1,'ljForce',false,...
                                   'diffusionForce',true,'ljPotentialWidth',0.1,'ljPotentialDepth',0.1);
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
        h.Step(chainPos); % move the histones                
        
        % Apply bending elasticity forces for beads inside the beam
        inBeam = r.handles.classes.domain.InDomain(chainPos,2);
        if any(inBeam)            
            bendingElasticityForce = 100 ;
            bForce = ForceManager.GetBendingElasticityForce(true,chainPos,r.objectManager.connectivity,bendingElasticityForce ,[]);
            % zero out forces outside the beam 
            bForce(~inBeam,:) = 0;
            % update the position of the chain 
            chainPos(inBeam,:) = chainPos(inBeam,:) + bForce(inBeam,:)*cp(1).dt;
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