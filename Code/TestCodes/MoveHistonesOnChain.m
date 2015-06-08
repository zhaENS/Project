function MoveHistonesOnChain
% This test function moves histones on a Rouse chain
close all
% % create chain and domain and register them in the ObjectManager

% -Relaxation time for the Rouse chain defined by the longest relaxation time--
% relaxation times - 300 beads ~= 1000 steps
%                    400 beads ~= 2000 steps
%                    500 beads ~= 3000 steps

% (numBeads*b)^2 / (3*D*pi^2) %[Doi &  Edwards p.96 eq. 4.37)
% using: b  = sqrt(3)~=1.7
%        dt = 0.01
%        D  = 1; diffusion const.
% (500*sqrt(3))^2 /(3*pi^2 * 1)
numRelaxationSteps = 1;
numRecordingSteps  = 200;
numBeamSteps       = 1000;
saveConfiguration  = false;
loadConfiguration  = false;

% Define the ROI for density estimation
rectX      = -4;
rectY      = -4;
rectWidth  = 2*abs(rectX);
rectHeight = 2*abs(rectY);

if loadConfiguration
    r = LoadConfiguration(loadConfiguration);
    r.params.simulator.showSimulation = true;
    r.InitializeGraphics
else
    % Initialize simulator framework parameters
    simulatorParams = SimulationFrameworkParams('showSimulation',true,...
                                                'numSteps',numRelaxationSteps,...
                                                'dimension', 3,...
                                                'dt',0.01,...
                                                'objectInteraction',false);                                                
                                                                                        
    % create an open domain
    openSpaceForces = ForceManagerParams('lennardJonesForce',false,...
                                         'LJPotentialWidth',0.1,...
                                         'LJPotentialDepth',0.1,...
                                         'diffusionForce',true,...
                                         'diffusionConst',1,...
                                         'mechanicalForce',false,...
                                         'mechanicalForceDirection','out',...
                                         'mechanicalForceCenter',[0 0 0],...
                                         'mechanicalForceMagnitude',0,...
                                         'dt',simulatorParams.simulator.dt);
    
    dp(1)         = DomainHandlerParams('domainShape','open',...
                                        'reflectionType','off',...
                                        'forceParams',openSpaceForces,...
                                        'domainWidth',100,...
                                        'dimension',simulatorParams.simulator.dimension);
                                                                                                                
    % % create a chain
    chainForces = ForceManagerParams('dt',simulatorParams.simulator.dt,...
                                     'springForce',false,...
                                     'bendingElasticityForce',false,...
                                     'bendingConst',0,...
                                     'springConst',1*simulatorParams.simulator.dimension*openSpaceForces.diffusionConst/(sqrt(3))^2,...
                                     'minParticleEqDistance',0);
    
    cp          = ChainParams('numBeads',300,...
                              'initializeInDomain',3,...
                              'forceParams',chainForces,...
                              'b',1*sqrt(3));
                          
    % create a cylindrical Beam as a domain
    cylinderForces = ForceManagerParams('diffusionForce',false,...
                                        'lennardJonesForce',false,...
                                        'morseForce',false);
    
    dp(2)          = DomainHandlerParams('domainShape','cylinder',...
                                         'reflectionType','off',...
                                         'domainWidth',2,...
                                         'domainHeight',70,...
                                         'forceParams',cylinderForces);
    

%     % create a sphere for visualization
    gSphereForce = ForceManagerParams('lennardJonesForce',false,...
                                      'diffusionForce',false,...
                                      'morseForce',false,...
                                      'mechanicalForce',false);
    
    dp(3)        = DomainHandlerParams('domainWidth',sqrt(cp.numBeads/6)*cp.b,...% radius of Gyration
                                       'domainCenter',[0 0 0],...
                                       'reflectionType','preserveEnergy',...
                                       'forceParams',gSphereForce);
    
    % register the object parameters in the simulator framework
    simulatorParams.SetDomainParams(dp);
    simulatorParams.SetChainParams(cp);
    
    % initialize simulator framework
    r = RouseSimulatorFramework(simulatorParams);
    % Run until relaxation time
    r.Run      
    % save configuration by name as chainPos
    SaveConfiguration(r,saveConfiguration);
end


% get the chain position to initialize the histone on it
[~,initialChainPosition] = r.objectManager.GetMembersPosition(1);
initialChainPosition     = initialChainPosition{1};

% Initialize histones with the chain position
histoneForce = ForceManagerParams('dt',r.params.simulator.dt,...
                                'diffusionForce',false,...
                                'diffusionConst',0,...
                                'mechanicalForce',false,...
                                'mechanicalForceDirection','out',...
                                'mechanicalForceMagnitude',0,...
                                'mechanicalForceCenter',[0 0 0],...
                                'lennardJonesForce',false,...
                                'LJPotentialWidth',0,...
                                'LJPotentialDepth',0);

histoneParams = HistoneParams('numHistones',1,'forceParams',histoneForce);
h             = Histone(histoneParams);

h.Initialize(initialChainPosition);

% Initialize graphics
if r.params.simulator.showSimulation
    [histHandle,pAxes,pHistHandle,pPolyHandle,~,projPlane2D,projPlane3D,dnaDensityHandle,histoneDensityHandle]= ...
        InitializeGraphics(rectX,rectY,rectWidth,rectHeight,r.simulationGraphics.handles.graphical.mainAxes, h.curPos,initialChainPosition);
end

sprintf('%s%f%s','Start recording at time: ',r.simulationData.step*r.params.simulator.dt,' sec.')

% Start recording densities with no beam effect
for sIdx = 1:numRecordingSteps    
       
    % advance one step
    [r,h,chainPos] = Step(r,h);
    
    % Update projectionPlane position according to the cm of the chain
    [rectX,rectY] = UpdateProjectionPlanePositionByCM(chainPos,rectWidth,rectHeight);
    
    % Calculate histone and DNA densities in the ROI
    [histoneDensity, dnaDensity] = CalculateDensitiesInROI(chainPos,h.curPos,rectX,rectY,rectWidth,rectHeight,histoneParams.numHistones);
    
    % update graphics
    if r.params.simulator.showSimulation
            UpdateGraphics(r.simulationData.step,r.params.simulator.dt,h.curPos,chainPos,histoneDensity,dnaDensity,...
            rectX,rectY,rectWidth, rectHeight,...
            dnaDensityHandle,histoneDensityHandle,histHandle,projPlane2D,projPlane3D,pHistHandle,pPolyHandle)
    end
end

% start beam
chainPos = r.objectManager.GetPosition(1);
chainPos = chainPos{1};
% move the beam to the chain's center of mass
UpdateBeamPosition(chainPos,r,2);
if r.params.simulator.showSimulation
    UpdateBeamGraphics(r,2)
end

% Create DNA damages in the ray
inBeam     = r.handles.classes.domain.InDomain(chainPos,2);
inBeamTemp = inBeam;
for iIdx = 2:numel(inBeam)-1
if inBeam(iIdx)
    inBeamTemp(iIdx+1) = true;
    inBeamTemp(iIdx-1) = true;
end
end
inBeam      = inBeamTemp;
inBeamInds  = find(inBeam);
% choose a fraction of inBeam to exclude
frInBeam      = randperm(sum(inBeam));
fracToExclude = 0; % fraction of damages to induce on the edges in the ray
inBeam(inBeamInds(frInBeam(1:round(numel(frInBeam)*fracToExclude))))= false;
affectedBeadsHandle = line('XData',chainPos(inBeam,1),'YData',chainPos(inBeam,2),...
                           'Parent',pAxes,'Marker','o','MarkerFaceColor','r','LineStyle','none');

r.runSimulation = true;

connectivityMat        = r.objectManager.GetConnectivityMapAsOne(1);
bendingElasticityConst = 1/(2*r.params.simulator.dt);

sprintf('%s%f%s','Beam shot at time: ',r.simulationData.step*r.params.simulator.dt,' sec.')

while all([r.simulationData.step<(r.simulationData.step+numBeamSteps),r.runSimulation])
    
    % Advance one simulation step
    [r,h,chainPos] = Step(r,h);
    [chainPos]     = ApplyDamageEffect(chainPos,inBeam,connectivityMat,bendingElasticityConst,r.params.simulator.dt);
    r.objectManager.DealCurrentPosition(1,chainPos)
    set(affectedBeadsHandle,'XData',chainPos(inBeam,1),'YData',chainPos(inBeam,2));
    % Update projectionPlane position according to the cm of the chain
    [rectX,rectY] = UpdateProjectionPlanePositionByCM(chainPos,rectWidth,rectHeight);
    
    % Calculate histone and DNA densities in the ROI
    [histoneDensity, dnaDensity] = CalculateDensitiesInROI(chainPos,h.curPos,rectX,rectY,rectWidth,rectHeight,histoneParams.numHistones);
    if r.params.simulator.showSimulation
%         gyrationsphereData   = r.simulationGraphics.handles.graphical.domain(3);% update gyration sphere center to the chain's center of mass
%         UpdateGyrationSphere(r,3,gyrationsphereData,chainPos);
        UpdateGraphics(r.simulationData.step,r.params.simulator.dt,h.curPos,chainPos,histoneDensity,dnaDensity,...
            rectX,rectY,rectWidth, rectHeight,dnaDensityHandle,histoneDensityHandle,histHandle,projPlane2D,projPlane3D,pHistHandle,pPolyHandle)
    end
end
end

function [r,h,chainPos] = Step(r,h)
r.Step;
[~,chainPos] = r.objectManager.GetMembersPosition(1);
chainPos     = chainPos{1};

% Move the histones
h.Step(chainPos,r.params.simulator.dt);% update current position
end

function [chainPos] = ApplyDamageEffect(chainPos,inBeam,connectivityMat,bendingElasticityConst,dt)
    % Apply forces on the chain falling in the beam
    bForce             = ForceManager.GetBendingElasticityForce(true,chainPos,connectivityMat,bendingElasticityConst ,[]);
    bForce(~inBeam,:)  = 0;
    % Update affected edges of the chain
    chainPos(inBeam,:) = chainPos(inBeam,:) + bForce(inBeam,:)*dt;
end

function [histHandle,pAxes,pHistHandle,pPolyHandle,dAxes,projPlane2D,projPlane3D,dnaDensityHandle,histoneDensityHandle]= InitializeGraphics(rectX,rectY,rectWidth,rectHeight,mAxes, histonePosition,initialChainPosition)
%     mAxes = r.simulationGraphics.handles.graphical.mainAxes;

% initialize histone graphics %TODO: incorporate histone graphics in the simulationGraphics class

histHandle = line('XData',histonePosition(:,1),...
                 'YData',histonePosition(:,2),...
                 'ZData',histonePosition(:,3),...
                 'marker','o',...
                 'MarkerFaceColor','y',...
                 'MarkerSize',10,...
                 'Parent',mAxes,...
                 'LineStyle','none');
daspect(mAxes,[1 1 1])
% create figure for the projection in the x-y plane
pFigure = figure;
pAxes   = axes('Parent',pFigure,...
    'XLim',get(mAxes,'XLim'), ...
    'YLim',get(mAxes,'YLim'),...
    'FontSize',30);
daspect(pAxes,[1 1 1])

% projected histone
pHistHandle  = line('XData',histonePosition(:,1),...
                    'YData',histonePosition(:,2),...
                    'Marker','o',...
                    'MarkerFaceColor','y',...
                    'MarkerSize',10,...
                    'Parent',pAxes,...
                    'LineStyle','none');
% projected Polymer
pPolyHandle  = line('XData',initialChainPosition(:,1),...
                    'YDAta',initialChainPosition(:,2),...
                    'Marker','o',...
                    'MarkerFaceColor','none',...
                    'MarkerEdgeColor','b',...
                    'markerSize',7,...
                    'Parent',pAxes,...
                    'LineStyle','-');

% create density axes
dFig   = figure;
dAxes  = axes('Parent',dFig,'NextPlot','add','Color','none','FontSize',30,'YLim',[0 1]);
xlabel(dAxes,'Time [sec]','FontSize',40);
ylabel(dAxes,'Density','FontSize',40)


% insert the projection to the 3d axes
projPlane3D = patch([rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],...
    [rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)],...
    'r', 'Parent',mAxes, 'FaceAlpha',0.5);
% insert the projection to the projection axes
projPlane2D = patch([rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)],...
    'r', 'Parent',pAxes, 'FaceAlpha',0.5);
dnaDensityHandle     = line('XData',0,'YData',NaN,'Parent',dAxes,'Color','b','LineWidth',4,'DisplayName','DNA density');
histoneDensityHandle = line('XData',0,'YData',NaN,'Parent',dAxes,'Color','y','Linewidth',4,'DisplayName','HistoneDensity');
legend(dAxes,get(dAxes,'Children'))


end

function [rectX,rectY]= UpdateProjectionPlanePositionByCM(chainPos,rectWidth, rectHeight)
% set x and y to be at the cm

rectX  = mean(chainPos(:,1));
rectY  = mean(chainPos(:,2));
rectX  = rectX- rectWidth/2;
rectY  = rectY- rectHeight/2;
end

function UpdateGraphics(step,dt,histPosition,chainPos,histoneDensity,dnaDensity,rectX,rectY,rectWidth, rectHeight,dnaDensityHandle,histoneDensityHandle,histHandle,projPlane2D,projPlane3D,pHistHandle,pPolyHandle)

UpdateDensityGraphics(step,dt,histoneDensity,dnaDensity,dnaDensityHandle,histoneDensityHandle)
UpdateHistoneGraphics(histHandle,histPosition)
UpdateProjectionPlaneGraphics(rectX,rectY, rectWidth, rectHeight,projPlane2D, projPlane3D)
UpdateProjectedHistoneGraphics(pHistHandle,histPosition);
UpdateProjectedPolymerGraphics(pPolyHandle,chainPos)
drawnow
end

function UpdateHistoneGraphics(histHandle,curPos)
set(histHandle,'XData',curPos(:,1),'YData',curPos(:,2),'ZData',curPos(:,3));
end

function UpdateProjectedHistoneGraphics(pHistHandle,curPos)
set(pHistHandle,'XData',curPos(:,1),'YData',curPos(:,2));
end

function UpdateProjectedPolymerGraphics(pPolyHandle,chainPos)
set(pPolyHandle,'XData',chainPos(:,1),'YData', chainPos(:,2))
end

function UpdateDensityGraphics(step,dt,histoneDensity,dnaDensity,dnaDensityHandle,histoneDensityHandle)
dnaDensityDataX     = get(dnaDensityHandle,'XData');
dnaDensityDataY     = get(dnaDensityHandle,'YData');
histoneDensityDataX = get(histoneDensityHandle,'XData');
histoneDensityDataY = get(histoneDensityHandle,'YData');

set(dnaDensityHandle,'XData',[dnaDensityDataX, step*dt], 'YData',[dnaDensityDataY,dnaDensity]);
% set(histoneDensityHandle,'XData',[histoneDensityDataX, step*dt], 'YData',[histoneDensityDataY,histoneDensity]);
%        line('XData',step*dt,'YData',histoneDensity,'Marker','o','markerFaceColor','y','Parent',dAxes)
%        line('XData',step*dt,'YData',dnaDensity,'Marker','o','Parent',dAxes,'markerFaceColor','b')
end

function UpdateProjectionPlaneGraphics(rectX,rectY, rectWidth, rectHeight,projPlane2D, projPlane3D)
% projPlane3d,projPlane2D are graphical handles
set(projPlane3D,'XData',[rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],...
    'YData',[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)]);
set(projPlane2D,'XData',[rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],...
    'YData',[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)]);
end

function UpdateGyrationSphere(simulationFrameworkHandle,domainNumber, gyrationSphereData,chainPos)
% update the sphere graphics
ccm  = mean(chainPos,1);% chain center of mass
xPos = gyrationSphereData.points.x;
yPos = gyrationSphereData.points.y;
zPos = gyrationSphereData.points.z;

xPos = xPos -xPos(1,1)+ccm(1);
yPos = yPos -yPos(1,1)+ccm(2);
zPos = zPos -zPos((size(zPos,1)-1)/2,1)+ccm(3);

set(gyrationSphereData.mesh,'XData',xPos,...
                            'YData',yPos,...
                            'ZData',zPos);
simulationFrameworkHandle.handles.classes.domain.params(3).domainCenter = ccm;
simulationFrameworkHandle.simulationGraphics.handles.graphical.domain(domainNumber).points.x = xPos;
simulationFrameworkHandle.simulationGraphics.handles.graphical.domain(domainNumber).points.y = yPos;
simulationFrameworkHandle.simulationGraphics.handles.graphical.domain(domainNumber).points.z = zPos;
end

function UpdateBeamPosition(chainPos,simulationFrameworkHandle,domainNumber)
% Change the parametr such that the inDomainfunction will work properly
cm = mean(chainPos,1);
simulationFrameworkHandle.handles.classes.domain.params(domainNumber).domainCenter = cm;
end

function UpdateBeamGraphics(simulationFrameworkHandle,domainNumber)
domainCenter = simulationFrameworkHandle.handles.classes.domain.params(domainNumber).domainCenter;
%     domainRad    = simulationFrameworkHandle.handles.classes.domain.params(domainNumber).domainWidth;
beamHandle   = simulationFrameworkHandle.simulationGraphics.handles.graphical.domain(domainNumber).mesh;
xPos         = get(beamHandle,'XData');
yPos         = get(beamHandle,'YData');
set(beamHandle,'XData',xPos+domainCenter(1),...
    'YData',yPos+domainCenter(2));

end

function [histoneDensity, dnaDensity] = CalculateDensitiesInROI(chainPos,histonePos,rectX,rectY, rectWidth, rectHeight,numHistones)
% calculate histone density
histoneDensity =  sum((histonePos(:,1)<=(rectX+rectWidth) & histonePos(:,1)>=rectX &...
    histonePos(:,2)<=(rectY+rectHeight) & histonePos(:,2)>=rectY))/numHistones;

% calculate DNA density
[dnaLengthIn,totalDNALength] = PolygonLengthInRoi(chainPos(:,1:2),rectX,rectY,rectWidth, rectHeight);
dnaDensity = dnaLengthIn./totalDNALength;

%         dnaDensity    = sum((chainPos(:,1)<=(rectX+rectWidth) & chainPos(:,1)>=rectX &...
%                              chainPos(:,2)<=(rectY+rectHeight) & chainPos(:,2)>=rectY));
end

function simFramework = LoadConfiguration(loadRelaxationConfiguration)

if loadRelaxationConfiguration
    [fName,fPath] = uigetfile('*.mat');
    load(fullfile(fPath,fName)); % simulation framework, saved as r
    f = whos;
    cInd = strcmpi({f.class},'RouseSimulatorFramework');
    if ~any(cInd)% class indicator
        error('this is not a file of class rouseSimulatorFramework')
    end
    simFramework = eval(f(cInd).name);
end
end

function SaveConfiguration(r,saveConfiguration)
% save RouseSimulatorFramework class
if saveConfiguration
    uisave('r')
end

end
