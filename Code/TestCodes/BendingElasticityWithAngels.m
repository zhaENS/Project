function force = BendingElasticityWithAngels(particlePosition, particleDistance,bendingConst,angle0)
%===========
% Calculate the force on each particle given its connectivity
numParticles = size(particlePosition,1);
dimension    = size(particlePosition,2);
force        = zeros(numParticles,dimension);

% cosTheta  = @(pos,dist,i) ((1./(dist(i,i+1)*dist(i+1,i+2)))*sum((pos(i,:)-pos(i+1,:)).*(pos(i+2,:)-pos(i+1,:))));

% F1        = @(pos,dist,i)(-1./dist(i+1,i+2))*((pos(i+1,:)-pos(i,:))./dist(i,i+1) + (1/dist(i+1,i+2)).*(pos(i+2,:)-pos(i+1,:)).*...
%                         cosTheta(pos,dist,i));
% F2        = @(pos,dist,i)(1/dist(i,i+1))*((pos(i+2,:)-pos(i+1,:))./dist(i+1,i+2) +...
%                 (pos(i+1,:)-pos(i,:)).*cosTheta(pos,dist,i)./dist(i,i+1));% d cos(theta_i)/dr_i
% F3        = @(pos,dist,i)  -(F1(pos,dist,i)-F2(pos,dist,i));

F1 = @(pos,dist,i,cosTheta0)((pos(i-2,:)-pos(i-1,:))*(pos(i,:)-pos(i-1,:))'./(dist(i-2,i-1)*dist(i-1,i)) -cosTheta0)*(pos(i-2,:)-pos(i-1,:))./(dist(i-2,i-1)*dist(i,i-1));
F2 = @(pos,dist,i,cosTheta0)((pos(i-1,:)-pos(i,:))*(pos(i+1,:)-pos(i,:))'./(dist(i-1,i)*dist(i+1,i)) -cosTheta0)* (2*pos(i,:)-pos(i-1,:)-pos(i+1,:))./(dist(i-1,i)*dist(i+1,i));
F3 = @(pos,dist,i,cosTheta0)((pos(i+2,:)-pos(i+1,:))*(pos(i,:)-pos(i+1,:))'./(dist(i+2,i+1)*dist(i+1,i)) -cosTheta0)* (pos(i+2,:)-pos(i+1,:))./(dist(i+2,i+1)*dist(i,i+1)); 

% % before 
% F1          = @(pos,dist,i)(pos(i-2,:)-pos(i-1,:))./(dist(i-2,i-1)*dist(i-1,i));
% % during 
% F2          = @(pos,dist,i)(2*pos(i,:)-pos(i-1,:)-pos(i+1,:))./(dist(i,i-1)*dist(i,i+1));
% % after 
% F3          = @(pos,dist,i)(pos(i+2,:)-pos(i+1,:))./(dist(i,i+1)*dist(i+1,i+2));
pos = particlePosition;
dist = particleDistance;

bAngle    = cos(angle0);

for pIdx = 1:numParticles    
%         inds = pIdx-2:pIdx;
%         inds = (inds>=1 & inds<= numParticles);
%         inds(2:3) = double(inds(3:-1:2));% replace 2 and 3;
%         
%         force(pIdx,:) = bendingConst*((cosTheta(particlePosition,particleDistance,pIdx-2)- bAngle)*F1(particlePosition,particleDistance,pIdx-2)*inds(1)+...
%                                       (cosTheta(particlePosition,particleDistance,pIdx)  - bAngle)*F2(particlePosition,particleDistance,pIdx)*inds(2)+...
%                                       (cosTheta(particlePosition,particleDistance,pIdx-1)- bAngle)*F3(particlePosition,particleDistance,pIdx-1)*inds(3));
%         
    if pIdx ==1                
%         force(pIdx,:) = bendingConst*(cosTheta(particlePosition,particleDistance,pIdx)-bAngle)*...
%                                    F3(particlePosition,particleDistance,pIdx);
%                                     F2(particlePosition,particleDistance,pIdx);
         force(pIdx,:) = bendingConst*F3(pos,dist,1,bAngle);
    elseif pIdx == 2 && (pIdx~=(numParticles-1));
        
%         force(pIdx,:) = bendingConst*((cosTheta(particlePosition,particleDistance,pIdx)-bAngle)*F3(particlePosition,particleDistance,pIdx)+...
%                                       (cosTheta(particlePosition,particleDistance,pIdx-1)-bAngle)*F2(particlePosition,particleDistance,pIdx));
         force(pIdx,:) = bendingConst*(F2(pos,dist,2,bAngle)+F3(pos,dist,2,bAngle));
        
    elseif pIdx == (numParticles-1) && (pIdx~=2)
%         force(pIdx,:) = bendingConst*((cosTheta(particlePosition,particleDistance,pIdx-2)-bAngle)*F1(particlePosition,particleDistance,pIdx)+...
%                                       (cosTheta(particlePosition,particleDistance,pIdx-1)-bAngle)*F2(particlePosition,particleDistance,pIdx));
          force(pIdx,:) = bendingConst*(F2(pos,dist,pIdx,bAngle)+F1(pos,dist,pIdx,bAngle));
    elseif pIdx == numParticles        
%         force(pIdx,:) = bendingConst*(cosTheta(particlePosition,particleDistance,pIdx-2)-bAngle)*F1(particlePosition,particleDistance,pIdx);        
         force(pIdx,:) = bendingConst*(F1(pos,dist,pIdx,bAngle));
    else
%         
%         force(pIdx,:) = bendingConst*((cosTheta(particlePosition,particleDistance,pIdx-2)- bAngle)*F1(particlePosition,particleDistance,pIdx)+...
%                                       (cosTheta(particlePosition,particleDistance,pIdx-1)- bAngle)*F2(particlePosition,particleDistance,pIdx)+...
%                                       (cosTheta(particlePosition,particleDistance,pIdx)  - bAngle)*F3(particlePosition,particleDistance,pIdx));
         force(pIdx,:)= bendingConst*(F1(pos,dist,pIdx,bAngle)+F2(pos,dist,pIdx,bAngle)+F3(pos,dist,pIdx,bAngle));
        
    end
    
end

force= -force;

end
