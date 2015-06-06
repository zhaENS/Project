% scr test bending elasticity with force at all joints
function scrTestBendingElasticityWithJointForces
close all

numPoints = 13;
dimension = 3;
diffusionConst = 1e-6;
rMat           =RouseMatrix(numPoints);
fixedParticles = 1;
% numSteps  = 3500;
alpha     = 1; % force parameter between 0 and 1
dt        = 0.1; % time step
% set the particles
particles = randn(numPoints,dimension);
if dimension ==2
    particles(:,3) = 0;
end

% set graphics
mainFig  = figure('Units','norm');
lim = [-1 1]*3.2;
mainAxes = axes('Parent',mainFig,'Color','k','Units','norm','XLim',lim, 'YLim',lim,'ZLim',lim);
uicontrol('Parent',mainFig,'Units','norm','Position',[0.05 0.05, 0.1 0.1],'String','stop','Callback',@StopButtonAction);
daspect(mainAxes,[1 1 1]);
setappdata(0,'stopButton',true)
cameratoolbar(mainFig);
l = line('XData',particles(:,1),'YData',particles(:,2),'ZData',particles(:,3),...
    'Marker','o','MarkerFaceColor','g','LineStyle','-','Parent',mainAxes,'Color','y');

while getappdata(0,'stopButton')
    
    forceVector = zeros(numPoints,3); % the force vector on each particle
    for pIdx = 2:numPoints-1
        midPoint = 0.5.*(particles(pIdx-1,:)+ particles(pIdx+1,:));
        
        % Calculate the force vectors on each particle
        
        % mid particle:
        forceVector(pIdx,:) = alpha.*(midPoint-particles(pIdx,:))+forceVector(pIdx,:);
        
        % neighboring particles
        
        % vector connecting neighbors
        vR = (particles(pIdx+1,:) - particles(pIdx,:)); % from left to right
        vL =  particles(pIdx-1,:) - particles(pIdx,:);
                
        % the normalized force vector of the neighbors
        dirR = (particles(pIdx+1,:) - midPoint);
        dirR = dirR./norm(dirR);
        
        dirL = (particles(pIdx-1,:) - midPoint);
        dirL = dirL./norm(dirL);
        
        a  = sqrt(sum((midPoint -particles(pIdx+1,:)).^2));
        b = norm(vR);
        if a<b
            t= (b-a)*alpha;
        else
            t = (b-a +b*alpha);
        end
        
        forceVector(pIdx+1,:) = forceVector(pIdx+1,:) + dirR*t ;
        
        b    = norm(vL);
        if a<b
            t= (b-a)*alpha;
        else
            t = (b-a +b*alpha);
        end
        forceVector(pIdx-1,:) = forceVector(pIdx-1,:) + dirR*t;
    end
    forceVector(fixedParticles,:) = 0;% keep first particle fixed
    springForce = -dimension*diffusionConst*rMat*particles;
    springForce(fixedParticles,:) = 0;
    diffusionForce = sqrt(2*diffusionConst*dt)*randn(size(particles));
    diffusionForce(fixedParticles,:) = 0;
    particles = particles+forceVector*dt +springForce*dt+ diffusionForce;
    set(l,'XData',particles(:,1),'YDAta',particles(:,2),'ZData',particles(:,3));
    drawnow
%     pause(dt)
end
end

function StopButtonAction(varargin)
if getappdata(0,'stopButton')
    setappdata(0,'stopButton', false)
else
    setappdata(0,'stopButton', true);
end
end