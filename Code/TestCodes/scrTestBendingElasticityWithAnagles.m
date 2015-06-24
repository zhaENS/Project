% scrTestBendingElasticityWithAnagles
% close all
numParticles     = 200;
dimension        = 2;

numSteps         = 1500;
dt               = 0.1;
angle0           = 1*pi;
b                = sqrt(3);
diffusionConst   = 1;
bendingConst     = 1*dimension*diffusionConst./b^2;
springConst      = 1*dimension*diffusionConst./b^2;
particlePosition = cumsum(sqrt(2*diffusionConst*dt)*randn(numParticles,3));
particlePosition(:,(1:3)>dimension) = 0;

connectivityMap    = (diag(ones(1,numParticles-1),1)+diag(ones(1,numParticles-1),-1))~=0;
% form a looped polymer 
% connectivityMap(1,end) = true;
% connectivityMap(end,1) = true;

minParticleDist    = b;
fixedParticleNum   = [];
afectedBeadsNumber = [];

% flags
diffusionFlag      = false;
springsFlag        = true;
bendingFlag        = true;

% rMat = RouseMatrix(numParticles);
gr           = sqrt(numParticles/6)*b^2 ;
% graphics
mainFig      = figure('Units','norm');
mainAxes     = axes('Parent',mainFig,'Units','norm','Color','k','XLim',[-gr gr],'YLim',[-gr gr],'ZLim',[-gr gr]);
dcmFigure    = figure('Units','norm'); 
dcmAxes      = axes('Parent',dcmFigure,'FontSize',25);xlabel(dcmAxes,'Time [sec]','FontSize',25); ylabel(dcmAxes,'mean dist. from cm','FontSize',25);
cameratoolbar(mainFig);
daspect(mainAxes,[1 1 1]);
particleHandle = line('XData',particlePosition(:,1),...
                      'Ydata',particlePosition(:,2),...
                      'Zdata',particlePosition(:,3),...
                      'Marker','o','MarkerFaceColor','g','MarkerEdgeColor','g',...
                      'LineStyle','-','Color','w','LineWidth',4,'Parent',mainAxes);
                  
% calculate the mean distance from cm 
cm            = mean(particlePosition,1);
dcm           = pdist2(particlePosition,cm);
meanDistCM    = line('XData',dt,'YData', mean(dcm),'Parent',dcmAxes,'Linewidth',3);
bLine         = line('XData',[0 0],'YData',[0 0],'color','g','Parent',dcmAxes);

affectedBeads = false(numParticles,1);
n = randperm(numParticles);
affectedBeads(n(1:10)) = true;
% particleDist  = b*ones(numParticles);
for sIdx = 1:numSteps
    
    particleDist = ForceManager.GetParticleDistance(particlePosition);
    
    % calculate bending elasticity force
    bendingForce = BendingElasticityWithAngels(particlePosition, particleDist,bendingConst,angle0);
    bendingForce(~affectedBeads,:) = 0;
    diffusionForce = sqrt(2*0.001*dt)*randn(numParticles,3);
    diffusionForce(:,(1:3)>dimension) =0;
%     springForce    = springConst*rMat*particlePosition;
    springForce   = ForceManager.GetSpringForce(springsFlag,particlePosition,particleDist,springConst,connectivityMap,minParticleDist,fixedParticleNum);
    
    particlePosition = particlePosition +bendingFlag*bendingForce*dt +springsFlag*springForce*dt+ diffusionFlag*diffusionForce;
                
    set(particleHandle,'XData', particlePosition(:,1),...
                       'YData',particlePosition(:,2),...
                       'ZData',particlePosition(:,3));
    
    cm  = mean(particlePosition,1);
    dcm = pdist2(particlePosition,cm);
    set(meanDistCM,'XData',[get(meanDistCM,'XData'), sIdx*dt],...
                   'YData',[get(meanDistCM,'YData'), mean(dcm)]);
    drawnow    
    
    if sIdx ==numSteps/2
                set(bLine,'XData',[sIdx sIdx].*dt,...
                          'YData',[0 max(get(meanDistCM,'YData'))],'Color','g','Parent',dcmAxes)
    end
    
    if sIdx>numSteps/2
        % turn off bending 
        bendingConst = bendingConst*(1 -10/(numSteps));%1*dimension*diffusionConst./(0.5*b)^2;
%         affectedBeads = false(numParticles,1);
        sprintf('%s','bending stoped') 
    end
end



