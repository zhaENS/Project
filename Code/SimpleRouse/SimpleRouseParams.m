classdef SimpleRouseParams<handle
% Parameters for the SimpleRouse class and simulations framework 
properties (SetObservable,GetObservable,AbortSet)
affineBeadsNum@double   % a fixed pair indices of affine beads
analyzeResults@logical  % perform analysis post simulations
b@double                % std of bead distance
beta@double             % ~~!!for beta~=2 the polymer is a beta polymer!!~~to set all betas to the same value insert one value only
calculateMSD@logical    % indicate whether to calculate the MSD for each bead (slows down simulations)
connectedBeads@double   % the beads connected other than the trivial connections. given by bead pairs
defaultRecipe@char      % default recipe filename 
diffusionConst@double   % difussion constant 
dimension@double        % dimension 1-3
dt@double               % time step for simulation
encounterDist@double    % bead encounter distance
kOff@double             % the detachment rate of affine beads
noiseCycle@double       % after how many steps we need to recalculate noise [unused]
noiseSTD@double         % std of noise added to steps of the simulation 
numBeads@double         % number of beads in the chain 
numRounds@double        % number of simulation rounds 
numSimulations@double   % number of simulations in each round
numSteps@double         % n times the number of steps until relaxation of the chain 
plot@logical            % TODO: should call the plotter  [obsolete]
recipeFileName@char     % recipe filename
recipeFolder@char       % recipe folder
recordPath@logical      % record bead position (slows down simulations)
saveBeadDist@char       % options [last/current/all/meanSquare] ( note that only for 'last' and 'all' the encounters can reliably be calculated)
saveResults@logical     % save results to folder (logical)  
springConst@double      % can be a scalar or a matrix the size of (numBeads) X (numBeads)
stiffConnectors@double  % fixed bead indices pairs for the number of stiff connectors 
end

methods
    
function obj = SimpleRouseParams
% class constructor
obj.SetDefaultParams;
% obj.CreateParamsListeners
end

function SetDefaultParams(obj)

% default params
% add listener to param change 
obj.numRounds       = 1;  % number of simulation rounds 
obj.numSimulations  = 20000; % number of simulations in each round
obj.dimension       = 3;
obj.numBeads        = 64;
obj.b               = sqrt(3);
obj.diffusionConst  = 1;
obj.connectedBeads  = []; % the beads connected other than the trivial connections. given by bead pairs
obj.encounterDist   = obj.b/5;
obj.dt              = 0.1*(obj.b^2)/(12*obj.diffusionConst);
obj.noiseSTD        = sqrt(2*obj.diffusionConst*obj.dt);
obj.noiseCycle      = 10000; % after how many steps we need to recalculate noise [unused]
obj.springConst     = -(obj.dimension*obj.diffusionConst*obj.dt/obj.b^2)*ones(obj.numBeads); % can be a scalar or a matrix the size of (numBeads) X (numBeads)

% set spring constant for the connected beads
for cIdx = 1:size(obj.connectedBeads,1)
    obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2))=obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2));
    obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1))=obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1));
end

d                   = sqrt(2*obj.diffusionConst*obj.dt);
obj.numSteps        = (.1)*(obj.b^2)/(6*(d^2)*sin(pi/(2*obj.numBeads))^2); % n times the number of steps until relaxation of the chain 
obj.numSteps        = round(obj.numSteps);
obj.stiffConnectors = [];    % fixed bead indices pairs for the number of loops
obj.affineBeadsNum  = [];    % a fixed pair indices of affine beads
obj.kOff            = 0.1;   % the detachment rate of affine beads
obj.beta            = 2;     % ~~!!for beta~=2 the polymer is a beta polymer!!~~to set all betas to the same value insert one value only
obj.recipeFolder    = fullfile(pwd,'SimpleRouse','Recipes');
obj.recipeFileName  = 'simpleRouseDebugRecipe';
obj.defaultRecipe   = 'simpleRouseDebugRecipe'; % default recipe file name
obj.saveBeadDist    = 'last'; % [last/current/all/meanSquare] ( note that only for 'last' and 'all' the encounters can reliably be calculated)
obj.calculateMSD    = true;  % indicate whether to calculate the MSD for each bead (slows down simulations)
obj.plot            = false;  % TODO: should call the plotter  [obsolete]
obj.recordPath      = false;  % record bead position (slows down simulations)
obj.analyzeResults  = true;  % perform analysis post simulations 
obj.saveResults     = true;  

end

end

methods (Access=private)
    
function CreateParamsListeners(obj)
    f= fieldnames(obj);
    for fIdx = 1:numel(f)
        addlistener(obj,f{fIdx},'PostSet',@obj.ParamChangeCallback);
    end
end

function ParamChangeCallback(obj,src,varargin)
        switch src.Name
            case 'b'
                obj.dt              = 0.1*(obj.b^2)/(12*obj.diffusionConst);
                obj.encounterDist   = obj.b/5;
            case 'dt'
                obj.noiseSTD        = sqrt(2*obj.diffusionConst*obj.dt);
                obj.numSteps        = round((0.02)*(obj.b^2)/(6*(sqrt(2*obj.diffusionConst*obj.dt)^2)*sin(pi/(2*obj.numBeads))^2)); % n times the number of steps until relaxation of the chain 
                obj.springConst     = -(obj.dimension*obj.diffusionConst*obj.dt/obj.b^2)*ones(obj.numBeads); % can be a scalar or a matrix the size of (numBeads) X (numBeads)
                % set spring constant for the connected beads
                for cIdx = 1:size(obj.connectedBeads,1)
                    obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2))=obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2));
                    obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1))=obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1));
                end
            case 'dimension'
                obj.springConst     = -(obj.dimension*obj.diffusionConst*obj.dt/obj.b^2)*ones(obj.numBeads); % can be a scalar or a matrix the size of (numBeads) X (numBeads)
                % set spring constant for the connected beads
                for cIdx = 1:size(obj.connectedBeads,1)
                    obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2))=obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2));
                    obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1))=obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1));
                end
            case 'numBeads' 
                obj.springConst     = -(obj.dimension*obj.diffusionConst*obj.dt/obj.b^2)*ones(obj.numBeads); % can be a scalar or a matrix the size of (numBeads) X (numBeads)
                % set spring constant for the connected beads
                for cIdx = 1:size(obj.connectedBeads,1)
                    obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2))=obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2));
                    obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1))=obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1));
                end
                obj.numSteps        = round((0.02)*(obj.b^2)/(6*(sqrt(2*obj.diffusionConst*obj.dt)^2)*sin(pi/(2*obj.numBeads))^2)); % n times the number of steps until relaxation of the chain 
%             case 'diffusionConst'
%                 obj.dt              = 0.1*(obj.b^2)/(12*obj.diffusionConst);                
%             case 'connectedBeads'
%                 obj.springConst     = -(obj.dimension*obj.diffusionConst*obj.dt/obj.b^2)*ones(obj.numBeads); % can be a scalar or a matrix the size of (numBeads) X (numBeads)
%                 % set spring constant for the connected beads
%                 for cIdx = 1:size(obj.connectedBeads,1)
%                     obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2))=obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2));
%                     obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1))=obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1));
%                 end
%             case 'springConst'
% %                 obj.springConst     = -(obj.dimension*obj.diffusionConst*obj.dt/obj.b^2)*ones(obj.numBeads); % can be a scalar or a matrix the size of (numBeads) X (numBeads)
%                 % set spring constant for the connected beads
%                 for cIdx = 1:size(obj.connectedBeads,1)
%                     obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2))=obj.springConst(obj.connectedBeads(cIdx,1), obj.connectedBeads(cIdx,2));
%                     obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1))=obj.springConst(obj.connectedBeads(cIdx,2), obj.connectedBeads(cIdx,1));
%                 end
        end
end

end
end

