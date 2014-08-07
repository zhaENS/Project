% scrAnalyzeTADDdata
% Analyze the data for TAD D alone.
close all
axesScale    = 'log'; % options [Linear\log];
analysisType = 'oneSided'; %options [oneSide\twoSided]% display the
if strcmpi(analysisType,'twoSided')
    axesScale = 'linear';
end

[fValues,data,title] = xlsread('D:\Ofir\Work\Copy\Polymer Chain Dynamics\Data\Luca\s65_contact_av.xlsx');
freq  = cell(size(fValues,1),1);
xData = freq;
for fIdx = 1:size(fValues,1)
    if strcmpi(analysisType,'oneSided')
        before = fValues(fIdx,fIdx:-1:1);
        after  = fValues(fIdx,fIdx+1:end);
        if numel(before)>numel(after)
            before(1:numel(after)) = before(1:numel(after))+after;
            freq{fIdx} = before;
        else
            after(1:numel(before)) = after(1:numel(before))+before;
            freq{fIdx} = after;
        end
        xData{fIdx} = 1:numel(freq{fIdx});
    elseif strcmpi(analysisType,'twoSided')
        freq{fIdx} = fValues(fIdx,:);
        xData{fIdx}      = -fIdx+1:size(fValues,2)-fIdx;
    end
end


a= axes('XScale',axesScale,'YScale',axesScale,'NextPlot','Add');
for fIdx = 1:numel(freq)
    line('XData',xData{fIdx},'YData',freq{fIdx},'Color',rand(1,3),'Parent',a);
end
title('display of normalized encounter frequency vs. distance [beads] of TAD D region');
xlabel(sprintf('%s%s',axesScale,' Distance in Beads'));
if strcmpi(axesScale,'log')
    ylabel('log(normalized # Encounter events)');
else
    ylabel('normalized #Encounter events');
end
line('XData',1:107,'YData',1./[1:107].^1.5,'Color','r','LineWidth',5)