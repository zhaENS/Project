function rouseParams        = SimpleRouseParams
% parameters for the SimpleRouse class and simulations
rouseParams.numRounds       = 5;% number of simulation rounds 
rouseParams.numSimulations  = 10;% number of simulations in each round
rouseParams.dimension       = 3;
rouseParams.numBeads        = 64;
rouseParams.b               = 0.1;
rouseParams.diffusionConst  = 1;
rouseParams.encounterDist   = rouseParams.b/10;
rouseParams.dt              = 0.2*(rouseParams.b^2)/(12*rouseParams.diffusionConst);
rouseParams.noiseSTD        = sqrt(2*rouseParams.diffusionConst*rouseParams.dt);
rouseParams.noiseCycle      = 10000;% after how many steps we need to recalculate noise [unused]
rouseParams.springConst     = -rouseParams.dimension*rouseParams.diffusionConst*rouseParams.dt/rouseParams.b^2;
d                           = sqrt(2*rouseParams.diffusionConst*rouseParams.dt);
numSteps                    = (0.5)*(rouseParams.b^2)/(6*(d^2)*sin(pi/(2*rouseParams.numBeads))^2); % n times the number of steps until relaxation of the chain 
rouseParams.numSteps        = round(numSteps);
rouseParams.stiffConnectors = []; % fixed bead indices pairs for the number of loops
rouseParams.connectedBeads  = []; % the beads connected other than the trivial connections. given by bead pairs
rouseParams.affineBeadsNum  = []; % a fixed pair indices of affine beads
rouseParams.kOff            = 0.1;% the detachment rate of affine beads
rouseParams.recipeFolder    = fullfile(pwd,'SimpleRouse','Recipes');
rouseParams.recipeFileName  = 'simpleRouseSimulateOneLoopWithInternalConnections';
rouseParams.defaultRecipe   = 'simpleRouseDebugRecipe'; % default recipe file name
rouseParams.saveBeadDist    = 'last'; % [last/current/all/meanSquare] ( note that only for 'last' and 'all' the encounters can reliably be calculated)
rouseParams.plot            = false;  % TODO: should call the plotter  [obsolete]
rouseParams.recordPath      = false;  % record bead position (slows down simulations)
rouseParams.analyzeResults  = true;  % perform analysis post simulations 
end