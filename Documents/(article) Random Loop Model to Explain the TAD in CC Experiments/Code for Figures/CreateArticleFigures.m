% this is a script to create the figures of the article 

% load the data for One TAD (this is a local folder) 
load('D:\Ofir\Work\ENS\PolymerChainDynamicsResults\simpleRouse_simpleRouseSimulateOneTADwithTail_ZeroTo30Loops_18_6_04-Oct-2014.mat')

% Figure 3
% positions of the boxes 
opx = [0.01 0.35 0.65 0.01 0.35 0.65 0.01 0.35 0.65]; 
opy = [0.72 0.72 0.72 0.39 0.39 0.39 0.05 0.05 0.05];

t= [1 5:5:30];
figure,
for tIdx = 1:numel(t)
 subplot(3,3,tIdx),
 imagesc(srf.encounterHistogram(:,:,t(tIdx))), 
 colormap hot, 
 daspect([1 1 1]),
 hold on , 
 text(250, 20, num2str(t(tIdx)),'Color','w','FontSize',20); % add number of loops 
 set(gca,'FontSize',25,'Position', [opx(tIdx) opy(tIdx) 0.27 0.27])
end

% load the data for two TADs
% Figure 4
load('D:\Ofir\Work\ENS\PolymerChainDynamicsResults\simpleRouse_simpleRouseSimulateTwoTADs0To30RandomLoops_18_11_02-Oct-2014.mat');
opx = [0.01 0.35 0.65 0.01 0.35 0.65 0.01 0.35 0.65]; 
opy = [0.72 0.72 0.72 0.39 0.39 0.39 0.05 0.05 0.05];
figure,
for tIdx = 1:numel(t)
 subplot(3,3,tIdx),
 imagesc(srf.encounterHistogram(:,:,t(tIdx))), 
 colormap hot, 
 daspect([1 1 1]),
 hold on , 
 text(250, 20, num2str(t(tIdx)),'Color','w','FontSize',20); % add number of loops 
 set(gca,'FontSize',25,'Position', [opx(tIdx) opy(tIdx) 0.27 0.27])
end

b= zeros(1,30);
for bIdx = 1:30
    b(bIdx) = srf.fitResults.mean{bIdx}.beta;
end
figure, 
plot(1:30,b,'-o','MarkerSize', 10,'MarkerFaceColor','g','LineWidth',7),
set(gca,'FontSize',50)