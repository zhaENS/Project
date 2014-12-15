function rouseParams        = SimpleRouseParams
% Parameters for the SimpleRouse class and simulations
rouseParams.numRounds       = 1; % number of simulation rounds 
rouseParams.numSimulations  = 4000; % number of simulations in each round
rouseParams.dimension       = 3;
rouseParams.numBeads        = 32;
rouseParams.b               = sqrt(3);
rouseParams.diffusionConst  = 1;
rouseParams.connectedBeads  = []; % the beads connected other than the trivial connections. given by bead pairs
rouseParams.encounterDist   = rouseParams.b/5;
rouseParams.dt              = 0.1*(rouseParams.b^2)/(12*rouseParams.diffusionConst);
rouseParams.noiseSTD        = sqrt(2*rouseParams.diffusionConst*rouseParams.dt);
rouseParams.noiseCycle      = 10000; % after how many steps we need to recalculate noise [unused]
rouseParams.springConst     = -(rouseParams.dimension*rouseParams.diffusionConst/rouseParams.b^2)*ones(rouseParams.numBeads); % can be a scalar or a matrix the size of (numBeads) X (numBeads)

% set spring constant for the connected beads
for cIdx = 1:size(rouseParams.connectedBeads,1)
    rouseParams.springConst(rouseParams.connectedBeads(cIdx,1), rouseParams.connectedBeads(cIdx,2))=rouseParams.springConst(rouseParams.connectedBeads(cIdx,1), rouseParams.connectedBeads(cIdx,2));
    rouseParams.springConst(rouseParams.connectedBeads(cIdx,2), rouseParams.connectedBeads(cIdx,1))=rouseParams.springConst(rouseParams.connectedBeads(cIdx,2), rouseParams.connectedBeads(cIdx,1));
end

d                           = sqrt(2*rouseParams.diffusionConst*rouseParams.dt);
numSteps                    = (0.01)*(rouseParams.b^2)/(6*(d^2)*sin(pi/(2*rouseParams.numBeads))^2); % n times the number of steps until relaxation of the chain 
rouseParams.numSteps        = round(numSteps);
rouseParams.stiffConnectors = [];    % fixed bead indices pairs for the number of loops
rouseParams.affineBeadsNum  = [];    % a fixed pair indices of affine beads
rouseParams.kOff            = 0.1;   % the detachment rate of affine beads
rouseParams.beta            = 2;     % ~~!!for beta~=2 the polymer is a beta polymer!!~~to set all betas to the same value insert one value only
rouseParams.recipeFolder    = fullfile(pwd,'SimpleRouse','Recipes');
rouseParams.recipeFileName  = 'simpleRouseCompareSpringConstantsBetweenChainAndRing';
rouseParams.defaultRecipe   = 'simpleRouseDebugRecipe'; % default recipe file name
rouseParams.saveBeadDist    = 'last'; % [last/current/all/meanSquare] ( note that only for 'last' and 'all' the encounters can reliably be calculated)
rouseParams.calculateMSD    = true;  % indicate whether to calculate the MSD for each bead (slows down simulations)
rouseParams.plot            = false;  % TODO: should call the plotter  [obsolete]
rouseParams.recordPath      = false;  % record bead position (slows down simulations)
rouseParams.analyzeResults  = true;  % perform analysis post simulations 
rouseParams.saveResults     = true;  


% change spring constant to matrix 

end