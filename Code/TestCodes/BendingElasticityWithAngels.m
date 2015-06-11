function force = BendingElasticityWithAngels(particlePosition, particleDistance,bendingConst,angle0)
%===========
% Calculate the force on each particle given its connectivity
numParticles = size(particlePosition,1);
dimension    = size(particlePosition,2);
force        = zeros(numParticles,dimension);

cosTheta  = @(pos,dist,i) ((1./(dist(i,i+1)*dist(i+1,i+2)))*sum((pos(i,:)-pos(i+1,:)).*(pos(i+2,:)-pos(i+1,:))));

% F1        = @(pos,dist,i)(-1./dist(i+1,i+2))*((pos(i+1,:)-pos(i,:))./dist(i,i+1) + (1/dist(i+1,i+2)).*(pos(i+2,:)-pos(i+1,:)).*...
%                         cosTheta(pos,dist,i));
% F2        = @(pos,dist,i)(1/dist(i,i+1))*((pos(i+2,:)-pos(i+1,:))./dist(i+1,i+2) +...
%                 (pos(i+1,:)-pos(i,:)).*cosTheta(pos,dist,i)./dist(i,i+1));% d cos(theta_i)/dr_i
% F3        = @(pos,dist,i)  -(F1(pos,dist,i)-F2(pos,dist,i));

% before 
F1          = @(pos,dist,i)(pos(i-2,:)-pos(i-1,:))./(dist(i-2,i-1)*dist(i-1,i));
% during 
F2          = @(pos,dist,i)(2*pos(i,:)-pos(i-1,:)-pos(i+1,:))./(dist(i,i-1)*dist(i,i+1));
% after 
F3          = @(pos,dist,i)(pos(i+2,:)-pos(i+1,:))./(dist(i,i+1)*dist(i+1,i+2));


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
        force(pIdx,:) = bendingConst*(cosTheta(particlePosition,particleDistance,pIdx)-bAngle)*...
                                   F3(particlePosition,particleDistance,pIdx);
%                                     F2(particlePosition,particleDistance,pIdx);
    elseif pIdx == 2 && (pIdx~=(numParticles-1));
        
        force(pIdx,:) = bendingConst*((cosTheta(particlePosition,particleDistance,pIdx)-bAngle)*F3(particlePosition,particleDistance,pIdx)+...
                                      (cosTheta(particlePosition,particleDistance,pIdx-1)-bAngle)*F2(particlePosition,particleDistance,pIdx));
        
    elseif pIdx == (numParticles-1)
        force(pIdx,:) = bendingConst*((cosTheta(particlePosition,particleDistance,pIdx-2)-bAngle)*F1(particlePosition,particleDistance,pIdx)+...
                                      (cosTheta(particlePosition,particleDistance,pIdx-1)-bAngle)*F2(particlePosition,particleDistance,pIdx));
        
    elseif pIdx == numParticles        
        force(pIdx,:) = bendingConst*(cosTheta(particlePosition,particleDistance,pIdx-2)-bAngle)*F1(particlePosition,particleDistance,pIdx);        
    else
        
        force(pIdx,:) = bendingConst*((cosTheta(particlePosition,particleDistance,pIdx-2)- bAngle)*F1(particlePosition,particleDistance,pIdx)+...
                                      (cosTheta(particlePosition,particleDistance,pIdx-1)- bAngle)*F2(particlePosition,particleDistance,pIdx)+...
                                      (cosTheta(particlePosition,particleDistance,pIdx)  - bAngle)*F3(particlePosition,particleDistance,pIdx));
        
    end
    
end

force= -force;

end
