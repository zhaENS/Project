% simulate a chain with affine bead. affine beads are defined in pairs in
% two separate regions of the chain

profile off

% Load model params
rouseParams = SimpleRouseParams;

% Preallocation
mEncounterProb = zeros(1,rouseParams.numBeads-1,rouseParams.numRounds);
encounterHist  = zeros(rouseParams.numBeads,rouseParams.numBeads,rouseParams.numRounds);
fitParams      = cell(1,rouseParams.numRounds);
gof            = cell(1,rouseParams.numRounds);
fittedExp      = zeros(1,rouseParams.numRounds);
fig  = figure('Units','norm');
hold on
ax       = axes('Parent',fig,'XScale','log','YScale','log');
fitModel = fittype('a*x.^(-b)');


for pIdx = 1:rouseParams.numRounds
    
    beadDist = zeros(rouseParams.numBeads,rouseParams.numBeads,rouseParams.numSimulations);
    
    % Create classs
    for sIdx = 1:rouseParams.numSimulations
       
        rouseChain = SimpleRouse(rouseParams);
        rouseChain.Initialize% run Rouse
        rouseChain.Run;
        beadDist(:,:,sIdx) = rouseChain.beadDist(:,:,end);
        sprintf('%s%s%s\n', 'simulation ', num2str(sIdx),' is done.')
    end
    
    % Sum over all experiments
    encounterHist(:,:,pIdx)  = sum(beadDist<rouseParams.encounterDist,3);% take the mean encounter for each bead pair
    encounterHist(:,:,pIdx) = encounterHist(:,:,pIdx)-diag(diag(encounterHist(:,:,pIdx)));
    mEncounterHist = zeros(rouseParams.numBeads,rouseParams.numBeads-1);
    for bIdx =1:size(encounterHist(:,:,pIdx),2)-1
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
                f =encounterHist(1,2:end,pIdx);% let it be zero
            end
        end
        mEncounterHist(bIdx,1:numel(f))=f;
    end
    
     % Calculate mean encounter probability  
    mEncounterProb(1,:,pIdx) = mean(mEncounterHist);
  
    % Plot the mean encounter probability and  the expected theoretical curve
    mh = max(mEncounterProb(1,:,pIdx));
    
    % Fit a line to the mean encounter probability from simulations
    [fitParams{pIdx}, gof{pIdx}] = fit((1:rouseParams.numBeads-1)',mEncounterProb(1,:,pIdx)',fitModel,...
        'StartPoint',[max(mEncounterProb(1,:,pIdx)),1.5],...
        'Robust','Bisquare');
    fittedExp(pIdx) = fitParams{pIdx}.b;
    % plot the encounter prob and fit 
     c = rand(1,3);
    line('XData',1:rouseParams.numBeads-1,...
        'YData',mEncounterProb(1,:,pIdx),...
        'Color',c,...
        'DisplayName',[num2str(numLoops), ' affine pairs'],...
        'Parent',ax), 
    line('XData',1:rouseParams.numBeads-1,...
         'YData',fitParams{pIdx}.a*((1:rouseParams.numBeads-1).^(-fitParams{pIdx}.b)),...
         'Color',c,...
         'DisplayName',['a=', num2str(fitParams{pIdx}.a),', b=', num2str(fitParams{pIdx}.b)],...
         'Parent',ax)
    drawnow
      
    % increase the encounter distance
    rouseParams.encounterDist = rouseParams.encounterDist+0.01;% reset the list
end

% plot the change in the fitted exponent value with increased loops
figure, plot(0:rouseParams.numRounds-1,fittedExp,'o-'), xlabel('number of loops'), ylabel('fitted exp');
% create a montage out of the encounter matrices
figure,
for pIdx = 1:rouseParams.numRounds
    subplot(5, 4, pIdx),imagesc(encounterHist(:,:,pIdx)), 
    colormap hot, 
    title([num2str(pIdx-1),' loops'])    

end


