% test eig fitting 
function SmoothUsingHeatEquation(sig)
% eMat is the encounter matrix 
% r       = RouseMatrix(size(eMat,1));
sig    = sig./sum(sig);% normalize the signal 
b0     = zeros(1,2*numel(sig));% initial guess
D      = 1;% diffusion constant 
lambda = (0:(2*numel(sig)-1))*pi./(2*numel(sig));
e1     = @(lambda,D,t)exp(-lambda.^2 *D*t);
[bs,fVal] = fminsearch(@(b)CoeffFinder(b,lambda,D,sig,e1),b0);

e2 = @(b,lambda,D,t)dot(b,exp(-lambda.^2 *D*t));
for tIdx=1:numel(sig)
solT(tIdx) = e2(bs,lambda,D,tIdx);
end
figure, plot(1:numel(sig),sig,'b',1:numel(sig), solT,'r')
end

function sol = CoeffFinder(b,lambda,D,sig,e)

for tIdx=1:numel(sig)
 solT(tIdx) = (dot(b,e(lambda,D,tIdx))-sig(tIdx))^2;
end
sol = sum(solT);
end