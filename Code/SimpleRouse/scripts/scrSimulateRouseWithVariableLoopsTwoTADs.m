% simulate a chain with affine bead. affine beads are defined in pairs in
% two separate regions of the chain

profile off
% model params
rouseParams = SimpleRouseParams;
tadSize     = 32;
numLoops    = 8; % 
% Preallocation
mEncounterProb  = zeros(1,rouseParams.numBeads-1,rouseParams.numRounds);
encounterHist   = zeros(rouseParams.numBeads,rouseParams.numBeads,rouseParams.numRounds);
fitParams       = cell(1,rouseParams.numRounds);
gof             = cell(1,rouseParams.numRounds);
fittedExp       = zeros(rouseParams.numBeads,rouseParams.numRounds);

figure,hold on,

for pIdx = 1:rouseParams.numRounds
    
    beadDist = zeros(rouseParams.numBeads,rouseParams.numBeads,rouseParams.numSimulations,rouseParams.numRounds);
    
    % Create classs
    for sIdx = 1:rouseParams.numSimulations
        
        % define connected beads for the first TAD
        r = randperm(tadSize);
        for rIdx= 1:numLoops
            rouseParams.connectedBeads(rIdx,1) = r(2*rIdx-1);
            rouseParams.connectedBeads(rIdx,2) = r(2*rIdx);
        end
        
       % define connected beads for the second TAD
        tadInds = rouseParams.numBeads:-1:rouseParams.numBeads-tadSize;
        r       = randperm(tadSize);
        for rIdx = 1:numLoops
            rouseParams.connectedBeads(numLoops+rIdx,1) = tadInds(r(2*rIdx-1));
            rouseParams.connectedBeads(numLoops+rIdx,2) = tadInds(r(2*rIdx));
        end
        
        % initialize Rosue class with the RouseParameters 
        rouseChain = SimpleRouse(rouseParams);
        rouseChain.Initialize
        rouseChain.Run% run Rouse
        beadDist(:,:,sIdx,pIdx) = rouseChain.beadDist(:,:,end);% save last distance matrix 
        sprintf('%s%s%s\n', 'simulation ', num2str(sIdx),' is done.')
    end
    
    % calculate encounter histogram 
    encounterHist(:,:,pIdx) = sum(beadDist<rouseParams.encounterDist,3);%  encounters 
    encounterHist(:,:,pIdx) = encounterHist(:,:,pIdx)-diag(diag(encounterHist(:,:,pIdx)));% remove main diagonal 
    
    %calculate  mean encounter histogram for each bead
    mEncounterHist = zeros(rouseParams.numBeads,rouseParams.numBeads-1);
    for bIdx =1:rouseParams.numBeads
            f = zeros(2,rouseParams.numBeads-1);
            f(1,1:numel(bIdx+1:rouseParams.numBeads)) = encounterHist(bIdx,bIdx+1:rouseParams.numBeads,pIdx);% right wing
            f(2,1:numel(1:bIdx-1)) = fliplr(encounterHist(bIdx,1:bIdx-1,pIdx));% flipped left wing
            f= sum(f);% sum left and eight 
            if sum(f) ~=0% normalize 
                normFac = ones(1,size(encounterHist(:,:,pIdx),2)-1);
                normFac(1:bIdx-1) = 2;
                f = f./normFac;
                f = f./sum(f);
            end
        mEncounterHist(bIdx,1:numel(f))=f;
    end
    

    % Calculate mean encounter probability  
    % the normalization term is determined by the number of terms avialable
    % for each bead. 
    normTerm = ones(1,rouseParams.numBeads-1);
    normTerm(1:floor(rouseParams.numBeads/2)) = rouseParams.numBeads;
    numTermsAtTail  = floor(rouseParams.numBeads/2)+1:rouseParams.numBeads-1;   
    normTerm(floor(rouseParams.numBeads/2)+1:end)= rouseParams.numBeads-2:-2:2;
    mEncounterProb(1,:,pIdx) = sum(mEncounterHist)./normTerm;
    
    % Plot the mean encounter probability and  the expected theoretical curve    
    % Fit a line to the Mean encounter probability from simulations
    hold on,
    fitModel = fittype('(1/sum(x.^-b)).*(x.^(-b))');
    fOptions = fitoptions(fitModel);
    fOptions.Robust = 'Bisquare';
    fOptions.Lower  = 0;
    fOptions.StartPoint = 1.5;
    dists               = (1:rouseParams.numBeads-1)';
    [fitParams{pIdx}, gof{pIdx}] = fit(dists,mEncounterProb(1,:,pIdx)',fitModel,...
        fOptions);            
    % plot the encounter prob and fit 
     c = rand(1,3);
    line('XData',dists,...
        'YData',mEncounterProb(1,:,pIdx),...
        'Color',c,...
        'DisplayName',[num2str(numLoops), ' affine pairs'],...
        'Parent',gca), 
    
    line('XData',1:rouseParams.numBeads-1,...
         'YData',(1/sum(dists.^(-fitParams{pIdx}.b))).*((dists).^(-fitParams{pIdx}.b)),...
         'Color',c,...
         'DisplayName',['a=', num2str(1/sum(dists.^(-fitParams{pIdx}.b))),', b=', num2str(fitParams{pIdx}.b)],...
         'Parent',gca)
      drawnow
   
  % fit the encounter function for each bead 
  for bdIdx = 1:rouseParams.numBeads
      beadFit = fit(dists,mEncounterHist(bdIdx,:)',fitModel,fOptions);
      fittedExp(bdIdx,pIdx) = beadFit.b;
  end
      
    % increase the number of affine pairs by 1 
%     numLoops = numLoops+1;
    rouseParams.connectedBeads = [];% reset the list
end

