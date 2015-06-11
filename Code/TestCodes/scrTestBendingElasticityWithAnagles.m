% scrTestBendingElasticityWithAnagles
close all
numParticles     = 10;
dimension        = 3;
particlePosition = cumsum(randn(numParticles,dimension));
numSteps       = 250;
dt             = 0.1;
angle0         = pi;
bendingConst   = -5.1;
diffusionConst = 0.1;
% flags
diffusionFlag      = false;
springsFlag        = false;
bendingFlag        = true;

b              = sqrt(3);
springConst  = -dimension*diffusionConst./b^2;
rMat = RouseMatrix(numParticles);
% graphics
mainFig      = figure('Units','norm');
mainAxes     = axes('Parent',mainFig,'Units','norm');
cameratoolbar(mainFig);
daspect(mainAxes,[1 1 1]);
particleHandle = line('XData',particlePosition(:,1),...
                      'Ydata',particlePosition(:,2),...
                      'Zdata',particlePosition(:,3),...
                      'Marker','o','MarkerFaceColor','k','LineStyle','-');
lastBeadHandle = line('XData',particlePosition(numParticles,1),...
                      'Ydata',particlePosition(numParticles,2),...
                      'Zdata',particlePosition(numParticles,3),...
                      'Marker','o','MarkerFaceColor','r',...
                      'MarkerEdgeColor','r','LineStyle','-','MarkerSize',7);

for sIdx = 1:numSteps
    
    particleDist = ForceManager.GetParticleDistance(particlePosition);
    % calculate bending elasticity force
    bendingForce = BendingElasticityWithAngels(particlePosition, particleDist,bendingConst,angle0);
    
    diffusionForce = sqrt(2*diffusionConst*dt)*randn(numParticles,dimension);
    springForce    = springConst*rMat*particlePosition;
    particlePosition = particlePosition +bendingFlag*bendingForce*dt +springsFlag*springForce*dt+ diffusionFlag*diffusionForce;
    
    set(particleHandle,'XData', particlePosition(:,1),...
        'YData',particlePosition(:,2),...
        'ZData',particlePosition(:,3));
     set(lastBeadHandle,'XData', particlePosition(numParticles,1),...
        'YData',particlePosition(numParticles,2),...
        'ZData',particlePosition(numParticles,3));
    drawnow
    pause(0.005)
end

