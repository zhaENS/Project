function scrCharacterizeMeanDensityInTADs
close all
% for TAD D
load('savedAnalysisTADD.mat');
e1j = a.encounterMatrix.average;
[m, meanB]=Analyze(e1j);
figure, plot(m),title(['TAD D, mean b= ', num2str(meanB)]), xlabel('bead'),ylabel('mean distance')

% TAD E
load('savedAnalysisTADE.mat')
e1j = a.encounterMatrix.average;
[m, meanB]=Analyze(e1j);
figure, plot(m), title(['TAD E, mean b=',num2str(meanB)]), xlabel('bead'),ylabel('mean distance')

% two TADS
load('savedAnalysisTADDAndE.mat')
e1j = a.encounterMatrix.average;
% e1j = e1j(80:250,80:250);
[m, meanB]=Analyze(e1j);
figure, plot(m), title(['TAD D and E. mean b=' num2str(meanB)]), xlabel('bead'),ylabel('mean distance')
% perform a moving averaging for the two TADS
width = 150;
m = zeros(1,size(e1j,1));
m(1:(width+1))= Analyze(e1j(1:(1+width), 1:(1+width)));% first step
for mIdx = 2:(size(e1j,1)-width)
[t,meanB] = Analyze(e1j(mIdx:(mIdx+width), mIdx:(mIdx+width)));
m(mIdx+width)= t(end);
temp = [m(mIdx:(mIdx+width));t];
for tIdx = 1:width
    f = ~isnan(temp(:,tIdx));
    meanVal = mean(temp(f,tIdx));
    if isempty(f)
        m(mIdx+tIdx-1)= NaN;
    else
        m(mIdx+tIdx-1)= meanVal;
    end
end
end
figure, plot(m), title(['moving average, TAD D and E. mean b=' num2str(meanB)]), xlabel('bead'),ylabel('mean distance')

end
function [m, meanB]=Analyze(e1j)
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
% for dIdx = 1:size(d,1); d(dIdx,d(dIdx,:)>(size(d,1)))=NaN; end

% estimate the mean density
m = zeros(1,size(e,1));
for dIdx = 1:size(d,1); m(dIdx)= mean(d(dIdx,~isnan(d(dIdx,:)))); end
% plot

end