function PiecewiseConstantSourceFourierSolution
% we display the series solution of the heat equation with a source 
% u_t=ku_xx+Q(x,t), with the initial conditions u(x,0) = const=c
% u_x(0,t) = u_x(L,t) = 0
% Q(x,t) is piece-wise constant ={1, t in [a_i,b_i], and 0 otherwise}
% the solution is evaluated at the point x=L/2 and hence becomes only a
% function of time 

n     = 150;% number of terms in the series
terms = cell(1,n);
for nIdx = 1:n
    nn = num2str(nIdx);
    lambda = sprintf('%s%s%s','(',num2str(nIdx),'*pi./L)');
    % the solution when the source is active, i.e Q=1
    eval(sprintf('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','termsQon{',nn,'}=@(k,t,L) exp(-',lambda,'*k*t).*((-1).^',nn,').*(2./(',nn,'*pi)).*(1-(-1).^',nn,').*((exp(',lambda,'*k*t)-1)./(',lambda,'*k)+1);'));
    % the solution when the source is inactive i.e Q=0
    eval(sprintf('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','termsQoff{',nn,'}=@(k,t,L) exp(-',lambda,'*k*t).*((-1).^',nn,').*(2./(',nn,'*pi)).*(1-(-1).^',nn,');'));
end
t  = 0:.01:15;
nt = numel(t);
k  = 0.1;
L  = 10;
% evaluate all the functions
solQon = zeros(n,nt);
for nIdx = 1:n
    solQon(nIdx,:) = termsQon{nIdx}(k,t,L);
    solQoff(nIdx,:) = termsQoff{nIdx}(k,t,L);
end
solTotalOn  = sum(solQon,1);
solTotalOff = sum(solQoff,1);
figure, plot(t,solTotalOn,'r',t,solTotalOff,'b')


end