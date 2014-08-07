% profile on 
% model params
numSimulations              = 1;
rouseParams.b               = 0.02;
rouseParams.numBeads        = 64;
rouseParams.dt              = rouseParams.b^3;
rouseParams.diffusionConst  = 1;
rouseParams.saveBeadDist    = 'current'; % [last/current/all/meanSquare] ( note that only for 'last' and 'all' the encounters can reliably be calculated)
rouseParams.plot            = false;
rouseParams.dimension       = 3;
rouseParams.recordPath      = true;
rouseParams.kOff            = 0.01;
rouseParams.encounterDist   = rouseParams.b/4;
rouseParams.stiffConnectors = [];
rouseParams.affineBeadsNum  = [];
rouseParams.connectedBeads  = [];
rouseParams.saveBeadPosition = false;

d        = sqrt(2*rouseParams.diffusionConst*rouseParams.dt);
numSteps =   2*(rouseParams.b^2)/(6*(d^2)*sin(pi/(2*rouseParams.numBeads))^2); % n times the number of steps until relaxation of the chain 
rouseParams.numSteps       = round(numSteps);


numAffinePairs = 10; 
% preallocation
beadDist = zeros(rouseParams.numBeads,rouseParams.numBeads,numSimulations);

% create classs
for sIdx = 1:numSimulations
    r = randperm(rouseParams.numBeads);
    rouseParams.affineBeadsNum(:,1) = r(1:numAffinePairs);
    rouseParams.affineBeadsNum(:,2) = r(numAffinePairs+1:2*numAffinePairs);
    rouseChain = SimpleRouse(rouseParams);
    rouseChain.Initialize% run Rouse 
    beadDist(:,:,sIdx) = rouseChain.beadDist(:,:,end);
     sprintf('%s%s%s\n', 'simulation ', num2str(sIdx),' is done.')
end

encounterHist = sum(beadDist<rouseParams.encounterDist,3);
mEncounterHist = zeros(rouseParams.numBeads,rouseParams.numBeads-1);

for bIdx =1:size(encounterHist,2)-1
    if bIdx~=1
        f = zeros(2,size(encounterHist,2)-1);
        f(1,1:numel(bIdx+1:size(encounterHist,2))) = encounterHist(bIdx,bIdx+1:size(encounterHist,2));
        f(2,1:numel(1:bIdx-1)) = fliplr(encounterHist(bIdx,1:bIdx-1));
        
        % divide by the number of elements 
        nz = zeros(1,size(f,2));
        nz(1,1:numel(1:bIdx-1))=1;
        nz(2,1:numel(bIdx+1:size(encounterHist,2)))=1;
        nz = sum(nz);
        f  = sum(f)./nz;
    else
        f= encounterHist(1,2:end);% the first line does not need averaging in both ends
    end
    mEncounterHist(bIdx,1:numel(f))=f;
end

% Calculate mean encounter probability 
mEncounterProb = zeros(1,size(mEncounterHist,2));
for mIdx = 1:size(mEncounterProb,2)
mEncounterProb(mIdx) = mean(mEncounterHist(~isnan(mEncounterHist(:,mIdx)),mIdx));
end

mEncounterProb = mEncounterProb./trapz(mEncounterProb);

figure, imagesc(encounterHist), colormap hot 
% h = encounterHist(1,2:end)./sum(encounterHist(1,2:end));

% plot the mean encounter probability and  the expected theoretical curve 
mh = max(mEncounterProb);
figure, plot(1:rouseParams.numBeads-1,mEncounterProb,'b',...
             1:rouseParams.numBeads-1,mh*(1:rouseParams.numBeads-1).^(-1.5),'r'),
xlabel('bead distance'), ylabel('encounter probability bead 1');

