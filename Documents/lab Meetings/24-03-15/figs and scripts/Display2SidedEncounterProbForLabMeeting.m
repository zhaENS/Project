% a script to create the graphs in the lab meeting presentation 
% lazzy loading 
if exist('a','var')&& isa(a,'AnalyzeEncounterFrequencies')     
     disp('data already loaded')        
else
    addpath(genpath(fullfile(pwd,'..','..','..','..','Code'))); % add the polymer chainDynamics folder
    % load the data
    load('ExperimentDataAnalysis\savedAnalysisTADDAndE.mat');    
end

%% figure 1
br.bead1=1:307; br.bead2=1:307; [eMat1p, eMat2p]= a.GetEncounterProbabilityMatrix(br,'average');
f= figure;
x=linspace(-306,306,size(eMat2p,2));
y=1:307;
surface(x,y,(eMat2p),'EdgeColor','none',...
    'FaceLighting','phong',...
    'FaceColor','interp',...
    'AmbientStrength',.3,...
    'SpecularExponent',2,...
    'BackFaceLighting','unlit',...
    'DiffuseStrength',.8);

currAxis= gca;

light('Parent',currAxis,...
    'Position',[-0.1 0.5 1],...
    'HandleVisibility','off');

camlight headlight, 
colormap('pink')
view(currAxis,[0 30]);
set(currAxis,'FontSize',30,'LineWidth',3,'XLim',[-306 306],'YLim',[1 307]);
xlabel('Distance from bead','FontSize',30);
ylabel('Bead num.','FontSize',30);
zlabel('encounter Prob.', 'Fontsize',30);
axis xy; grid on


%% figure 2

% choose at random numSig and display them 
r     = randperm(307);
figure,
numSig = 5;
for i=1:numSig
    subplot(numSig,1,i), plot(eMat2p(r(i),308:end),'color',rand(1,3),'LineWidth',3),
    title(['Bead ', num2str(r(i))],'FontSize',25)
    yl = get(gca,'YLim');
    set(gca,'FontSize',25,'LineWidth',2,'XLim',[1 150],'YLim',[-0.0 yl(2)])
end


