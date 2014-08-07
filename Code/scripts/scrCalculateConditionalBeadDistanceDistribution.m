% scrTestConditionalBeadDistanceDistribution
close all
minDist = 0.2; % contact distance
dists   = cell(64,63);
figure, hold on
for b1Idx = 1
    for b2Idx = b1Idx+1:64
        res       = results.beadDistMat(b1Idx,b2Idx,:);
        f         = find(res(:)<minDist);
        d         = results.beadDistMat(b1Idx,:,f);
        meanDist  = cell(1,numel(f));
        mds       = [];
        for dIdx = 1:numel(f)
            midLoop = ((b2Idx+b1Idx)/2);
            ma = max([1:b1Idx-1,b2Idx-ceil(midLoop-1),numel(b2Idx+1:64)]);
            mat1 = zeros(4,ma);
            mat1(1,1:numel(b1Idx-1:-1:1)) = d(1,b1Idx-1:-1:1);
            mat1(2,1:numel(b1Idx+1:ceil(midLoop-1)))=d(1,b1Idx+1:ceil(midLoop-1),dIdx);
            mat1(3,1:numel(b2Idx-1:-1:floor(midLoop+1))) = d(1,b2Idx-1:-1:floor(midLoop+1),dIdx);
            mat1(4,1:numel(b2Idx+1:64))=d(1,b2Idx+1:end,dIdx);
            numNonZeros = sum(mat1~=0);
            meanDist{dIdx} = sum(mat1)./numNonZeros;
            
        end
        
        for mIdx = 1:numel(f)
            mds(mIdx,:) = meanDist{mIdx};
        end
        
        dists{b1Idx,b2Idx} = mean(mds);
        line('Parent',gca,...
            'XData',1:numel(dists{b1Idx,b2Idx}),...
            'YData',dists{b1Idx,b2Idx},...
            'Color',rand(1,3),...
            'DisplayName',['bead',num2str(b1Idx),'LoopBead',num2str(b2Idx)])
        
        mds = [];
        A = [];
        B = [];
        C = [];
        D = [];
    end
end
% plot the expected curve
line('Parent',gca,...
    'XData',1:63,...
    'YData',0.2*(1:63).^0.5,...
    'Color','r',...
    'LineWidth',5)
