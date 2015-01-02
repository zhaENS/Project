function [cDelta,icDelta, detCdelta, v,lambda]= RouseRingMatrixAnalysis(numIncrements)
% Construct the covarience matrix of the increment of the Brownian bridge
% (rouse ring) 
N       = numIncrements;       % number of increments in the bridge 
cDelta  = zeros(N); % covariance matrix 
icDelta = zeros(N); % inverse covariance matrix
for s = 1:N
    for t= 1:N
        cDelta(s,t) = min([t+1,s+1])-min([t+1,s])-min([t,s+1])+min([s,t])-(2*N-1)/N^2;
    end
end

w= @(s,t,N) exp(2*pi*1i*t./N).^s;

% The eigenvectors 
v = zeros(N);
for vIdx = 1:N
    v(:,vIdx) = (1/sqrt(N))*w(vIdx-1,(0:N-1),N);
end

% The covariance is a circulant matrix, with eigenvalues 
lambda = zeros(N,1);
for nIdx = 1:N
    lambda(nIdx) = dot(cDelta(1,:),w(nIdx-1,(0:N-1),N));
end

% and determinant 
detCdelta = prod(lambda);

% Calculating the inverse covariance matrix using the eigenvalues and
% eigenvectors

for s = 1:N
    for t = 1:N
        icDelta(s,t) = (N/(1-N) + sum(w(s-1,(1:N-1),N).*w(t-1,(1:N-1),N)))/N;
    end
end


% show the pdf for each bead 
% x = zeros(1,N);
% x(1,1) = 1;
f=@(x,N,detCdelta,icDelta)  (1./((2*pi)^numel(x) * detCdelta))^(1/2) *exp(-0.5*(x*icDelta*x'));
figure, hold on
for i=1:10
 x = ones(1,N)*i;
plot(i,real(f(x,N,detCdelta,icDelta)),'.')
end



