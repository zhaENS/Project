% scrSimulate transiant loops
% at each simulation two beads are defined to have an affinity. the beads
% can connect if they are close and detach at some probability.
% at the end, the mean encounter matrix is calculated

numAffineBeads        = 25;
tadSize               = 50;
params.numBeads       = 128;
params.b              = 0.1;
params.diffusionConst = 1;
params.dt             = 1e-4;
params.dimension      = 3;
params.noiseSTD       = sqrt(2*params.diffusionConst*params.dt);
params.numSteps       = 30000;
params.affineBeadsNum = []; % two column vector of beads with affinity to one another
params.kOff           = 0.05; % [some units]
params.plot           = false;
params.encounterDist  = params.b/2;

step = 1;
for b1Idx = 1:numAffineBeads-2
    for b2Idx = b1Idx+2:numAffineBeads
        % define 2 pairs of affine beads
        params.affineBeadsNum = [b1Idx b2Idx; ...
            b1Idx+params.numBeads-tadSize-1,b2Idx+params.numBeads-tadSize];
        s = SimpleRouse(params);
        s.Initialize;
        results.beadsEncounterHist(:,:,step) = s.encounterHist;
        results.affineBeads{step} = params.affineBeadsNum;
        sprintf('%s%s%s%s%s%s%s%s%s%s\n','Simulation ', num2str(step),' done. ' ,...
            'Affine Beads TAD1:[', num2str(params.affineBeadsNum(1,:)),...
            ']. Affine Beads TAD2: [', num2str(params.affineBeadsNum(2,:)),...
            ']. SimulatioTime: ' , num2str(s.simulationTime),' [sec]')
        step = step+1;
    end
end

% calculate the mean encounter matrix over all simulations 
m = mean(results.beadsEncounterHist,3);
mbeaddist = zeros(params.numBeads);
for bIdx =1:size(m,2)-1
    if bIdx~=1
        f = zeros(2,size(m,2)-1);
        f(1,1:numel(bIdx+1:size(m,2))) = m(bIdx,bIdx+1:size(m,2));
        f(2,1:numel(1:bIdx-1)) = fliplr(m(bIdx,1:bIdx-1));
        nz = sum(f~=0);
        f  = sum(f)./nz;
    else
        f= m(1,2:end);
    end
    mbeaddist(bIdx,1:numel(f))=f;
end
