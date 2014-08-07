% scrSimulateInitial3Cconditions
startme
% profile on 
% model params
numSimulations             = 1000;
rouseParams.b              = 0.01;
rouseParams.numBeads       = 64;
rouseParams.dt             = rouseParams.b^3;
rouseParams.diffusionConst = 1;
rouseParams.saveBeadDist   = 'current'; % [last/current/all/meanSquare] ( note that only for last and all the encounters can reliably be calculated)
rouseParams.plot           = false;
rouseParams.dimension      = 3;
rouseParams.encounterDist  = rouseParams.b/2;
d        = sqrt(2*rouseParams.diffusionConst*rouseParams.dt);
numSteps =   2*(rouseParams.b^2)/(6*(d^2)*sin(pi/(2*rouseParams.numBeads))^2); % n times the number ofsteps until relaxation of the chain 
rouseParams.numSteps       = round(numSteps);

% preallocation
beadDist = zeros(rouseParams.numBeads,rouseParams.numBeads,numSimulations);

% create classs

for sIdx = 1:numSimulations
    rouseChain = SimpleRouse(rouseParams);
    rouseChain.Initialize% run rouse 
    beadDist(:,:,sIdx) = rouseChain.beadDist;
     sprintf('%s%s%s\n', 'simulation ', num2str(sIdx),' is done.')
end
encounterHist = sum(beadDist<rouseParams.encounterDist*3,3);
mEncounterHist = zeros(rouseParams.numBeads);
for bIdx =1:size(encounterHist,2)-1
    if bIdx~=1
        f = zeros(2,size(encounterHist,2)-1);
        f(1,1:numel(bIdx+1:size(encounterHist,2))) = encounterHist(bIdx,bIdx+1:size(encounterHist,2));
        f(2,1:numel(1:bIdx-1)) = fliplr(encounterHist(bIdx,1:bIdx-1));
        nz = sum(f~=0);
        f  = sum(f)./nz;
    else
        f= encounterHist(1,2:end);
    end
    mEncounterHist(bIdx,1:numel(f))=f;
end

mHist = mean(beadDist,3);
encounterHist = mHist<rouseParams.encounterDist;
encounterHist = encounterHist-eye(size(encounterHist)).*encounterHist;% remove the diagonal
figure, imagesc(encounterHist), colormap hot 
h = mEncounterHist(1,:)./sum(mEncounterHist(1,:));
mh = max(h);
figure, plot(1:63,h,'b',1:63,mh*(1:63).^(-1.5),'r'), xlabel('bead distance'), ylabel('encounter probability bead 1')
