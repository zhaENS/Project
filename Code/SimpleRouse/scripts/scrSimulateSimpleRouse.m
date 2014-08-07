% scrSimulateSimpleRouse
% In this script we run a simple Rouse chain with no constraints. The
% statistical properties of the chain are examined
% simulation 

profile off 
format long
% model params
rouseParams = SimpleRouseParams;
% preallocation 
encounterHist  = zeros(rouseParams.numBeads,rouseParams.numBeads,rouseParams.numRounds);
mEncounterProb = zeros(1,rouseParams.numBeads-1,rouseParams.numRounds);
for pIdx = 1:rouseParams.numRounds
% preallocation
beadDist = zeros(rouseParams.numBeads,rouseParams.numBeads,rouseParams.numSimulations);

% create classs
bd = cell(1,rouseParams.numSimulations);
for sIdx = 1:rouseParams.numSimulations
    rouseChain = SimpleRouse(rouseParams);
    rouseChain.Initialize
    rouseChain.Run;% run rouse 
    beadDist(:,:,sIdx) = rouseChain.beadDist(:,:,end);% take the last recorded bead distances
    sprintf('%s%s%s\n', 'simulation ', num2str(sIdx),' is done.')
end

% 
% % compare the mean square difference of bond length to the theoretical
% % result
% g   = @(D,dt,b) 54*(D*dt/b)^2 +4*D*dt;
% g(rouseParams.diffusionConst,rouseParams.dt,rouseParams.b)

% sum over all experiments
encounterHist(:,:,pIdx)  = sum(beadDist<rouseParams.encounterDist,3);          % take the mean encounter for each bead pair
encounterHist(:,:,pIdx)  = encounterHist(:,:,pIdx)-diag(diag(encounterHist(:,:,pIdx)));            % remove diagonal values
mEncounterHist = zeros(rouseParams.numBeads,rouseParams.numBeads-1); % preallocation

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

figure, imagesc(encounterHist(:,:,pIdx)), colormap hot 

% plot the mean encounter probability and  the expected theoretical curve 
mh = max(mEncounterProb(1,:,pIdx));
figure, plot(1:rouseParams.numBeads-1,mEncounterProb(1,:,pIdx),'b',...
             1:rouseParams.numBeads-1,mh*(1:rouseParams.numBeads-1).^(-1.5),'r'),
xlabel('bead distance'), ylabel('encounter probability bead 1');

% Fit a line to the mean encounter probability from simulations 
hold on, 
fitModel = fittype('1./(sum(x.^-b))*(x.^(-b))');
fitOptions = fitoptions(fitModel);
set(fitOptions,'Robust','Bisquare',...
                        'Lower',0,...
                        'StartPoint',1.5,...
                        'TolX',1e-8);
[fitParams, gof] = fit((1:rouseParams.numBeads-1)',mEncounterProb(1,:,pIdx)',fitModel,...
                        fitOptions);
dists = 1:rouseParams.numBeads-1;
plot(dists,(1/sum(dists.^-fitParams.b))*((dists).^(-fitParams.b)),'g')
end
