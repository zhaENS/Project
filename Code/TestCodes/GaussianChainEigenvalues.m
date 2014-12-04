function GaussianChainEigenvalues % unfinished
% This function follows the Eichinger algorithm to calculate the
% eigenvalues of generalized Gaussian chain.
% Sub units of the chain are identical 

numBeads   = 10; % number of beads
chainGraph = zeros(numBeads); % the connected compoent 
% label nodes 

[nodeLabels] = LabelNodes(chainGraph);


function [nodeLabels] = LabelGraph(chainGraph)
% Label the graph nodes, the output nodeLabels has the same size as the
% input chainGraph with represents the connectivity matrix of the graph 

% preallocation
nodeLabels = zeros(size(chainGraph));
% 1. find vertices with degree>2
s      = sum(chainGraph,2);
d      = find(s>2);
newLab = 1:numel(d);
diagC  = diag(chainGraph);
diagC(s>2) = newLab;
