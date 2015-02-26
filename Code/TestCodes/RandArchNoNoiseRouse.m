function RandArchNoNoiseRouse
% random architecture and no noise Rouse chain solution

r = RouseMatrix(10);
% r(1,1)  = 3;  r(1,10)  = -1; r(1,5)  = -1;
% r(5,1) = -1; r(5,5) = 4; r(5,7) = -1; 
% r(5,7) = -1; r(5,8) = -1; r(5,9) = -1; r(5,5) = 5;
% r(7,5) = -1; r(7,7) = 3;
% r(8,5) = -1; r(8,8) = 3;
% r(9,5) = -1; r(9,9) = 3;
r(1,10) = -1; r(1,1) = 2;
r(10,1) = -1; r(10,10) = 2;
[v,e] = eig(r);

p0 = (1:10)';% initial position
z0 = v*p0;   % in normal coordinates
k  = .4;     % spring constant
b  = 1;      % bead distance STD
t  = 0:20;   % time 
% solution of the equation with no noise 
for tIdx = 1:numel(t)
z(:,tIdx)   = z0.*exp(-(k/b^2)*diag(e)*t(tIdx)); %normal coordinates
pos(:,tIdx) = v'*z(:,tIdx);  % bead position 
end
f= figure;
a= axes('Parent',f);
for pIdx = 1:size(r,1)
 line('XData',t,'YData',pos(pIdx,:),'displayname',['bead ',num2str(pIdx)],'parent',a,'Color',rand(1,3))
end
hold on 
% plot(repmat(t(end),size(r,1),1),p0,'o')
legend(get(a,'Children'));
diag(e)