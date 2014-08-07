% In This script we model a chain with affine beads. the pair of affine
% beads goes from 0 to 16, filling the whole TAD.

% Load model params
rouseParams = SimpleRouseParams;
fitModel = fittype('a.*(x.^(-b))');
% preallocation
encounterHist = zeros(rouseParams.numBeads,rouseParams.numBeads,rouseParams.numRounds);
mEncounterProb = zeros(1,rouseParams.numBeads-1,rouseParams.numRounds);
numAffinePairs = 0; % start with no affine pair, (A control chain)

for pIdx = 1:rouseParams.numRounds
% preallocation
beadDist = zeros(rouseParams.numBeads,rouseParams.numBeads,rouseParams.numSimulations);

% create classs
for sIdx = 1:rouseParams.numSimulations
    if numAffinePairs~=0
     r = randperm(rouseParams.numBeads);    
     rouseParams.affineBeadsNum(:,1) = r(1:numAffinePairs);
     rouseParams.affineBeadsNum(:,2) = r(numAffinePairs+1:2*numAffinePairs);
    end
    rouseChain = SimpleRouse(rouseParams);
    rouseChain.Initialize% 
    rouseChain.Run;% run Simulation
    beadDist(:,:,sIdx) = rouseChain.beadDist(:,:,end);
     sprintf('%s%s%s\n', 'simulation ', num2str(sIdx),' is done.')
end

% sum over all experiments
encounterHist(:,:,pIdx)  = sum(beadDist<rouseParams.encounterDist,3);% take the mean encounter for each bead pair
mEncounterHist = zeros(rouseParams.numBeads,rouseParams.numBeads-1);
% Calculate the encounter probability over all beads
for bIdx =1:size(encounterHist,2)-1
        if bIdx~=1
            f = zeros(2,size(encounterHist,2)-1);
            f(1,1:numel(bIdx+1:size(encounterHist,2))) = encounterHist(bIdx,bIdx+1:size(encounterHist,2),pIdx);
            f(2,1:numel(1:bIdx-1)) = fliplr(encounterHist(bIdx,1:bIdx-1,pIdx));
            f= sum(f);
            if sum(f) ~=0
            f = f./sum(f);
            end

        else
            if sum(encounterHist(1,2:end,pIdx))~=0
                f= encounterHist(1,2:end,pIdx)./sum(encounterHist(1,2:end,pIdx));
            else
                f =encounterHist(1,2:end,pIdx);
            end
        end
        mEncounterHist(bIdx,1:numel(f))=f;
       
end


% Calculate mean encounter probability  
mEncounterProb(1,:,pIdx) = mean(mEncounterHist);

% plot the mean encounter probability and  the expected theoretical curve 
mh = max(mEncounterProb(1,:,pIdx));
line('XData',1:rouseParams.numBeads-1,...
     'YData',mEncounterProb(1,:,pIdx),...
     'Color','b',...
     'Parent',gca);
 
line('XData',1:rouseParams.numBeads-1,...
     'YData',mh*(1:rouseParams.numBeads-1).^(-1.5),...
     'Color','r',...
     'Parent',gca)


% Fit a line to the mean encounter probability from simulations 
[fitParams{pIdx}, gof] = fit((1:rouseParams.numBeads-1)',mEncounterProb(1,:,pIdx)',fitModel,...
                        'StartPoint',[mh,2/3],...
                        'Robust','Bisquare');

line('XData',1:rouseParams.numBeads-1,...
     'YData',fitParams{pIdx}.a.*((1:rouseParams.numBeads-1).^(-fitParams{pIdx}.b)),...
     'Color','g',...
     'Parent',gca)

 numAffinePairs = numAffinePairs+1;% increase the number of affine pairs
 rouseParams.affineBeadsNum = [];
end

xlabel('bead distance'), ylabel('encounter probability bead 1');