% scrTestBendingElasticityWithAnagles
close all
numParticles     = 25;
dimension        = 3;

numSteps         = 550;
dt               = 0.01;
angle0           = pi;
b                = sqrt(3);
diffusionConst   = 1;
bendingConst     = 10*dimension*diffusionConst./b^2;
springConst      = 5*dimension*diffusionConst./b^2;
particlePosition = cumsum(sqrt(2*diffusionConst*dt)*randn(numParticles,dimension));

connectivityMap  = (diag(ones(1,numParticles-1),1)+diag(ones(1,numParticles-1),-1))~=0;
minParticleDist  = 0;
fixedParticleNum = [];
afectedBeadsNumber  = [];

% flags
diffusionFlag      = false;
springsFlag        = true;
bendingFlag        = true;

% rMat = RouseMatrix(numParticles);
gr = 2*sqrt(numParticles/6)*b^2 ;
% graphics
mainFig      = figure('Units','norm');
mainAxes     = axes('Parent',mainFig,'Units','norm','Color','k','XLim',[-gr gr],'YLim',[-gr gr],'ZLim',[-gr gr]);
cameratoolbar(mainFig);
daspect(mainAxes,[1 1 1]);
particleHandle = line('XData',particlePosition(:,1),...
                      'Ydata',particlePosition(:,2),...
                      'Zdata',particlePosition(:,3),...
                      'Marker','o','MarkerFaceColor','g','MarkerEdgeColor','g',...
                      'LineStyle','-','Color','w','LineWidth',4);
lastBeadHandle = line('XData',particlePosition(numParticles,1),...
                      'Ydata',particlePosition(numParticles,2),...
                      'Zdata',particlePosition(numParticles,3),...
                      'Marker','o','MarkerFaceColor','r',...
                      'MarkerEdgeColor','r','LineStyle','-','MarkerSize',9);
                  

affectedBeadsHandle = line('XData',particlePosition(afectedBeadsNumber ,1),...
                           'Ydata',particlePosition(afectedBeadsNumber ,2),...
                           'Zdata',particlePosition(afectedBeadsNumber ,3),...
                          'Marker','o','MarkerFaceColor','g','MarkerEdgeColor','g',...
                           'LineStyle','-','Color','w','LineWidth',4);
 affectedBeads          = false(numParticles,1);
for sIdx = 1:numSteps
    
    particleDist = ForceManager.GetParticleDistance(particlePosition);
    % calculate bending elasticity force
    bendingForce = BendingElasticityWithAngels(particlePosition, particleDist,bendingConst,angle0);
    bendingForce(~affectedBeads,:) = 0;
    diffusionForce = sqrt(2*diffusionConst*dt)*randn(numParticles,dimension);
%     springForce    = springConst*rMat*particlePosition;
    springForce   = ForceManager.GetSpringForce(springsFlag,particlePosition,particleDist,springConst,connectivityMap,minParticleDist,fixedParticleNum);
    
    particlePosition = particlePosition +bendingFlag*bendingForce*dt +springsFlag*springForce*dt+ diffusionFlag*diffusionForce;
                
    set(particleHandle,'XData', particlePosition(:,1),...
        'YData',particlePosition(:,2),...
        'ZData',particlePosition(:,3));
     set(lastBeadHandle,'XData', particlePosition(numParticles,1),...
        'YData',particlePosition(numParticles,2),...
        'ZData',particlePosition(numParticles,3));
    
     set(affectedBeadsHandle,'XData', particlePosition(affectedBeads,1),...
        'YData',particlePosition(affectedBeads,2),...
        'ZData',particlePosition(affectedBeads,3));
    
    drawnow    
    
    if sIdx ==200
        affectedBeads          = false(numParticles,1);
        affectedBeads(afectedBeadsNumber ) = true;
        set(affectedBeadsHandle,'MarkerFaceColor','r','MarkerEdgeColor','r')
        
    end
end

