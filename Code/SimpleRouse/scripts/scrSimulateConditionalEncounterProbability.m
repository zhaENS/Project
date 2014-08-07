% Calculate the probability that bead A meets B before it meets C
% bead A is set to 1, bead C is set to 32, bead B is changed from 2 to 31.
profile off
format long
% model params
rouseParams = SimpleRouseParams;
% preallocation
encounterHist  = zeros(rouseParams.numBeads,rouseParams.numBeads,rouseParams.numRounds);
mEncounterProb = zeros(1,rouseParams.numBeads-1,rouseParams.numRounds);

beadA   = 1;
beadB   = 2;
beadC   = rouseParams.numBeads;
eFreqAB = zeros(1,rouseParams.numRounds);
eFreqAC = zeros(1,rouseParams.numRounds);
for pIdx = 1:rouseParams.numRounds    
    for sIdx = 1:rouseParams.numSimulations
        rouseChain = SimpleRouse(rouseParams);
        rouseChain.Initialize % run rouse
        keepRunning    = true;
        simulationStep = 1;
        while keepRunning
            % get the distances throughout the simulation
            rouseChain.Step % advance one simulation step
            dBeadAB     = rouseChain.beadDist(beadA,beadB,end);
            dBeadAC     = rouseChain.beadDist(beadA,beadC,end);
            encounterAB = dBeadAB(:)<rouseParams.encounterDist;
            encounterAC = dBeadAC(:)<rouseParams.encounterDist;
            if encounterAB&&~encounterAC
                eFreqAB(pIdx) = eFreqAB(pIdx)+ 1;
                keepRunning = false;
            elseif ~encounterAB&&encounterAC
                eFreqAC(pIdx) = eFreqAC(pIdx)+1;
                keepRunning = false;
            end
            if simulationStep>=rouseParams.numSteps
                keepRunning = false;
            else 
                simulationStep = simulationStep+1;
            end
        end
        sprintf('%s%s%s\n', 'simulation ', num2str(sIdx),' is done.')
    end
    % transform to probability
    eFreqAB(pIdx) = eFreqAB(pIdx)/rouseParams.numSimulations;
    eFreqAC(pIdx) = eFreqAC(pIdx)/rouseParams.numSimulations;
    
    %  advance the bead B by one
    beadB          = beadB+1;
end

figure, 
axes('NextPlot','add');

    line('XData',2:1+rouseParams.numRounds,...
             'YData',eFreqAB,...
             'Color','b',...
             'Parent',gca,...
             'DisplayName','prob A-B before A-C');
         
   line('XData',2:1+rouseParams.numRounds,...
             'YData',eFreqAC,...
             'Color','r',...
             'Parent',gca,...
             'DisplayName','prob A-C before A-B');      
  xlabel('bead B')
  ylabel('conditional Encounter probability')
  legend(get(gca,'Children')) 
  
         