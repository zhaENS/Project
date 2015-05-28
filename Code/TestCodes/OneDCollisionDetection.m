function [partCollide]=OneDCollisionDetection
close all 
% Define a polygon in 3D 
numPart  = 20;
numEdges = 20;
dt       = 0.1;
chain    = cumsum(randn(numEdges,3));

% parametrize the polygon according to arc-length (cumulative edge length)
edgeLength = sqrt(sum(diff(chain,1).^2,2));
p          = [0; cumsum(edgeLength)]; 

% divide by the total length 
p = p./p(end);

% choose random position on the chain 
prevPos = rand(numPart,1);

% move the particles on the 1D line 

curPos = prevPos+2*randc(numPart,1)*dt;

% plot 
mainFig  = figure;
mainAxes = axes('Parent',mainFig);
% plot parametize chain 
line('XData',p,'YData', zeros(1,numEdges));
% green - prevPos 
% red   - curPos 

% plot particle curPos
line('XData',curPos,'YData',zeros(1,numPart),'Marker','o','MarkerFaceColor','r','Parent',mainAxes)
text(curPos,0.1*ones(1,numPart),{1:numPart}','Parent',mainAxes,'Color','g')

% plot particles prevPos 
line('XData',prevPos,'YData',zeros(1,numPart),'Marker','o','MarkerFaceColor','g','Parent',mainAxes)
text(prevPos,0.1*ones(1,numPart),{1:numPart}','Parent',mainAxes,'Color','r')

% Collision detection 
% calculate the velocity of each particle 
partVelocity = (curPos-prevPos)./dt;
velocityDist = pdist2(partVelocity,partVelocity);
partDist     = pdist2(prevPos,prevPos);
t            = (partDist./velocityDist)./dt;

% look for t in the range [0 1]
collisionMat = (t>=0 & t<=1);
[partCollide(:,1), partCollide(:,2)] = find(triu(collisionMat));

% Reflect 

% transform back to 3D




end