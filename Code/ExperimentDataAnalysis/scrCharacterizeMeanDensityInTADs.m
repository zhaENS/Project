function scrCharacterizeMeanDensityInTADs
close all
% for TAD D
load('savedAnalysisTADD.mat');
e1j = a.encounterMatrix.average;
Analyze(e1j);
% TAD E
load('savedAnalysisTADE.mat')
e1j = a.encounterMatrix.average;
Analyze(e1j)

% two TADS
load('savedAnalysisTADDAndE.mat')
e1j = a.encounterMatrix.average;
Analyze(e1j)

end
function Analyze(e1j)
% estimate b from the encounter probability of nearest neighbors
e = zeros(size(e1j));
for dIdx = 1:size(e1j,1);e(dIdx,:) = e1j(dIdx,:)./sum(e1j(dIdx,~isnan(e1j(dIdx,:))));end
% get super and sub diagonal
b1 = e(diag(true(1,size(e,1)-1),1) | diag(true(1,size(e,1)-1),-1));
% remove zeros
b1= b1(b1>0);
% estimate the lognormal probability
[vars] = mle(b1,'distribution','logn');
meanB = exp(vars(1)+0.5*vars(2)^2);
% plug meanB to the distance estimate
d = ((3./(2*pi*(meanB^2))).^1.5)./e1j ;
d(isinf(d))=NaN;
% remove point bigger than the allowed distance
for dIdx = 1:size(d,1); d(dIdx,d(dIdx,:)>(size(d,1)-dIdx))=NaN; end
% estimate the mean density
m = zeros(1,size(e,1));
for dIdx = 1:size(d,1); m(dIdx)= mean(d(dIdx,~isnan(d(dIdx,:)))); end
% plot
figure, plot(m), xlabel('bead'),ylabel('mean distance')

end