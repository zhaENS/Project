% test eig fitting 
function SmoothUsingHeatEquation(sig)
% eMat is the encounter matrix 
% r       = RouseMatrix(size(eMat,1));
sig       = sig./sum(sig);% normalize the signal 
b0        = zeros(1,2*numel(sig));% initial guess
b0        = ones(1,30);
D         = 1;% diffusion constant 
% lambda    = (0:(2*numel(sig)-1))*pi./(2*numel(sig));
[bs,fVal] = fminsearch(@(b)SumOfExp(b,D,sig),b0);
% % e1        = @(lambda,D,t)exp(-lambda.^2 *D*t);
% [bs,fVal] = fminsearch(@(b)CoeffFinder(b,lambda,D,sig,e1),b0);
% 
e2 = @(b,lambda,lags,D,t)dot(b,exp(-lambda.^2 *D.*(lags-t)));
% 
solT = zeros(1,numel(sig));
for tIdx=1:numel(sig)
solT(tIdx) = e2(bs(1:10),bs(11:20),bs(21:30),D,tIdx);
end
figure, plot(1:numel(sig),sig,'b',1:numel(sig), solT,'r')
end


function sol= SumOfExp(b,D,sig)
s = zeros(1,numel(sig));
for tIdx =1:numel(sig)
s(tIdx) = (dot(b(1:10),exp(-b(11:20).^2 *D.*(b(21:30)-tIdx))-sig(tIdx)).^2);
end

sol = sum(s);
end

function sol = CoeffFinder(b,lambda,D,sig,e)

for tIdx=1:numel(sig)
 solT(tIdx) = (dot(b,e(lambda,D,tIdx))-sig(tIdx))^2;
end
sol = sum(solT);
end