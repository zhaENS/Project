% simulate a chain with affine bead. affine beads are defined in pairs in
% two separate regions of the chain

profile off
% model params
numRounds                   = 13;
numSimulations              = 1;
numLoops                    = 0;
tadSize                     = 27;
rouseParams.dimension       = 3;
rouseParams.b               = 0.1;
rouseParams.numBeads        = 64;
rouseParams.dt              = 1e-4;% changes later to prevent simulation 'explosion'
rouseParams.diffusionConst  = 1;
rouseParams.saveBeadDist    = 'last'; % [last/current/all/meanSquare] ( note that only for 'last' and 'all' the encounters can reliably be calculated)
rouseParams.plot            = false;

rouseParams.recordPath      = false;
rouseParams.encounterDist   = rouseParams.b/2;
rouseParams.stiffConnectors = [];
rouseParams.affineBeadsNum  = [];
rouseParams.kOff            = [];
rouseParams.connectedBeads  = [];

% choose k smaller than 1 
k = 0.005;
% determine dt such that the spring constant equals k
rouseParams.dt = k*(rouseParams.b^2)/rouseParams.diffusionConst;

d        = sqrt(2*rouseParams.diffusionConst*rouseParams.dt);
numSteps =   1.5*(rouseParams.b^2)/(6*(d^2)*sin(pi/(2*rouseParams.numBeads))^2); % n times the number of steps until relaxation of the chain
rouseParams.numSteps       = round(numSteps);



% Preallocation
mEncounterProb = zeros(1,rouseParams.numBeads-1,numRounds);
encounterHist  = zeros(rouseParams.numBeads,rouseParams.numBeads,numRounds);
fitParams      = cell(1,numRounds);
gof            = cell(1,numRounds);

figure,hold on,

for pIdx = 1:numRounds
    
    beadDist = zeros(rouseParams.numBeads,rouseParams.numBeads,numSimulations);
    
    % Create classs
    for sIdx = 1:numSimulations
        % define affine beads for the first TAD
        r = randperm(tadSize);
        for rIdx= 1:numLoops
            rouseParams.connectedBeads(rIdx,1) = r(2*rIdx-1);
            rouseParams.connectedBeads(rIdx,2) = r(2*rIdx);
        end
%         % define affine beads for the second TAD
%         tadInds = rouseParams.numBeads:-1:rouseParams.numBeads-tadSize;
%         r = randperm(tadSize);
%         for rIdx = 1:numLoops
%             rouseParams.connectedBeads(numLoops+rIdx,1) = tadInds(r(2*rIdx-1));
%             rouseParams.connectedBeads(numLoops+rIdx,2) = tadInds(r(2*rIdx));
%         end
        
        rouseChain = SimpleRouse(rouseParams);
        rouseChain.Initialize% run Rouse
        beadDist(:,:,sIdx) = rouseChain.beadDist(:,:,end);
        sprintf('%s%s%s\n', 'simulation ', num2str(sIdx),' is done.')
    end
    
    % Sum over all experiments
    encounterHist(:,:,pIdx)  = sum(beadDist<rouseParams.encounterDist,3);% take the mean encounter for each bead pair
    mEncounterHist = zeros(rouseParams.numBeads,rouseParams.numBeads-1);
    for bIdx =1:size(encounterHist(:,:,pIdx),2)-1
        if bIdx~=1
            f = zeros(2,size(encounterHist(:,:,pIdx),2)-1);
            f(1,1:numel(bIdx+1:size(encounterHist(:,:,pIdx),2))) = encounterHist(bIdx,bIdx+1:size(encounterHist,2),pIdx);
            f(2,1:numel(1:bIdx-1)) = fliplr(encounterHist(bIdx,1:bIdx-1,pIdx));
            
            % Divide by the number of elements
            nz = zeros(1,size(f,2));
            nz(1,1:numel(1:bIdx-1))=1;
            nz(2,1:numel(bIdx+1:size(encounterHist(:,:,pIdx),2)))=1;
            nz = sum(nz);
            f  = sum(f)./nz;
        else
            f= encounterHist(1,2:end,pIdx);
        end
        mEncounterHist(bIdx,1:numel(f))=f;
    end
    
    % Calculate mean encounter probability    
    for mIdx = 1:size(mEncounterHist,2)
        mEncounterProb(1,mIdx,pIdx) = mean(mEncounterHist(~isnan(mEncounterHist(:,mIdx)),mIdx));
    end
    
    mEncounterProb(1,:,pIdx) = mEncounterProb(1,:,pIdx)./trapz(1:rouseParams.numBeads-1,mEncounterProb(1,:,pIdx));
        
    % Plot the mean encounter probability and  the expected theoretical curve
    mh = max(mEncounterProb(1,:,pIdx));
    
    % Fit a line to the mean encounter probability from simulations
    hold on,
    fitModel = fittype('a.*(x.^(-b))');
    [fitParams{pIdx}, gof{pIdx}] = fit((1:rouseParams.numBeads-1)',mEncounterProb(1,:,pIdx)',fitModel,...
        'StartPoint',[max(mEncounterProb(1,:,pIdx)),1.5]);
    
    % plot the encounter prob and fit 
     c = rand(1,3);
    line('XData',1:rouseParams.numBeads-1,...
        'YData',mEncounterProb(1,:,pIdx),...
        'Color',c,...
        'DisplayName',[num2str(numLoops), ' affine pairs'],...
        'Parent',gca), 
    line('XData',1:rouseParams.numBeads-1,...
         'YData',fitParams{pIdx}.a.*((1:rouseParams.numBeads-1).^(-fitParams{pIdx}.b)),...
         'Color',c,...
         'DisplayName',['a=', num2str(fitParams{pIdx}.a),', b=', num2str(fitParams{pIdx}.b)],...
         'Parent',gca)
      drawnow
      
    % increase the number of affine pairs by 1 
    numLoops = numLoops+1;
   
end
