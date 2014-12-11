classdef AnalyzeEncounterFrequencies<handle
    % Experimental setting:
    % investigation of TAD E and D in Nora et al 2012.
    % the DNA segment consists of 920432 bp, which are cut using restiction
    % enzymes into restiction fragments, using HINDIII restriction enzyme.
    % there are forward and backward restriction segments. 124 forward
    % (FOR) and 126 reverse (REV) segments of variable size. HindIII
    % recognized the site AAGCTT and cleaves it between the AA, leaving a
    % sticky end
    
    properties
        allData
        beadData
        segmentData
        encounterMatrix
        classes
        peaks
        results
        params        
        %         fitModel     = fittype(@(slope,bias,x)(bias.*x.^(-slope)))
        creationDate = date;% class creation date      
    end
    
    events
    end
    
    methods
        function obj = AnalyzeEncounterFrequencies()
            % Class constructor
        end
        
        function Initialize(obj,params)
            % Initialization sequence             
            obj.LoadDefaultParams
            if exist('params','var')
                obj.SetInputParams(params)
            end
            obj.params.numBeads    = numel(unique([obj.params.beadRangeToAnalyze.bead1,obj.params.beadRangeToAnalyze.bead2]));
            obj.classes.peakFinder = PeakCalling;
            obj.ReadData;% load data from xls            
            obj.AnalyzeData
                                  
        end
        
        function LoadDefaultParams(obj)
            % loads default parameters
            % default location of xls file
            obj.params.xlsFilePath        = fullfile(pwd,'..','Data','Luca','MC_TAD_DE_E14d0_rep2+3_true.xlsx');
            obj.params.fillGapsBy         = 'sameValuesAsBoundary'; % when the data is  missing (i.e, the segment takes the position of few beads,
            obj.params.beadRangeToAnalyze.bead1 = (1:307); % if empty, analyze all,
            obj.params.beadRangeToAnalyze.bead2 = (1:307); % if empty, analyze all,
            obj.params.replicateName      = {'Rep1','Rep2','Average'};
            obj.params.beadSizeInbp       = 3000;      
            obj.params.fdrThresh          = 0.01; % threshhol for false discovery rate
            % fit options 
            
            obj.params.fitModel             = fittype(@(slope,x)(1./(sum(x.^(-slope)))).*x.^(-slope));
            obj.params.fOptions             = fitoptions(obj.params.fitModel);
            obj.params.fOptions.TolX        = 1e-8;
            obj.params.fOptions.TolFun      = 1e-8;
            obj.params.fOptions.MaxFunEvals = 1e3;
            obj.params.fOptions.MaxIter     = 1e3;
            obj.params.fOptions.StartPoint  = 1.5; % [slope]
            obj.params.fOptions.Lower       = 0.05; % [slope]
            obj.params.fOptions.Robust      = 'Bisquare';
            
            % statistics optimization params 
            obj.params.stOptions  = statset('Robust','on');
        end
        
        function SetInputParams(obj,params)
            % load input params. the param input must be a structure with
            % field names identical to the properties in obj.params
            fNames = fieldnames(params);
            for fIdx = 1:numel(fNames)
                if strcmpi(fNames{fIdx},'beadRangeToAnalyze')
                    obj.CheckAnalysisRangeMonotonicity(params.beadRangeToAnalyze)
                end
                
                obj.params.(fNames{fIdx}) = params.(fNames{fIdx});
            end
        end        
        
        function ReadData(obj)
            % a sequence of calls to list segments lengths and beads and
            % the transformations (interms of indices) between segments and
            % beads
            obj.ProcessAllData;
            
            obj.ProcessDataByBead;
            
            obj.ProcessDataBySegment;
            
        end
        
        function AnalyzeData(obj)
            
            obj.CreateEncounterFrequencyMatrices;
            fNames = lower(obj.params.replicateName);
            for fIdx = 1:numel(fNames)
             [oneSide, twoSides] = obj.ProcessEncounters(obj.params.beadRangeToAnalyze,fNames{fIdx});
              obj.beadData.encounterData.twoSides.(fNames{fIdx})= twoSides;
              obj.beadData.encounterData.oneSide.(fNames{fIdx}) = oneSide;
            end
            obj.FitData;
            
            obj.FitMeanModel;
            
            % Find peaks in the encounter data 
            obj.PeakCalling
        end
        
        function ProcessAllData(obj)
            % Read data from xls
            [values,~, obj.allData.all]= xlsread(obj.params.xlsFilePath);
            titles = {'beadStart1','beadEnd1','beadStart2','beadEnd2',...
                      'rep1Count','rep2Count','averageCount','stdCount'};
%             
%             % Read data only in the valid range of beads
%             if ~isempty(obj.params.beadRangeToAnalyze)
%                 range = obj.params.beadRangeToAnalyze;
%                 % take only beads in the range and such that have interactions with beads
%                 % in the range
%                 values = values(values(:,2)<=range(2) & values(:,1)>=range(1) & values(:,4)<=range(2) & values(:,3)>=range(1),:);
%             end
            
            % store the data for each replicate of experiment and the
            % average
            for tIdx = 1:numel(titles);
                obj.allData.(titles{tIdx})= values(:,tIdx);
            end
        end
        
        function ProcessDataByBead(obj)
            % List unique segments and segment encounter
            % The first index is the number of bead in which the restriction segment
            % starts. The second index (column) is the bead in which the
            % segment ends  
            allBeads1                     = unique([obj.allData.beadStart1,obj.allData.beadEnd1] ,'rows','stable');% unique segment 1
            allBeads2                     = unique([obj.allData.beadStart2,obj.allData.beadEnd2],'rows','stable'); % unique segment 2
            allBeads                      = unique([allBeads1;allBeads2],'rows');% uniqe segment encounter pairs
            obj.beadData.allInteractions  = allBeads;
                        
            obj.beadData.numBeads         = 307;% min([obj.params.beadRangeToAnalyze(2),307])- max([obj.params.beadRangeToAnalyze(1),1])+1;
            beadInRange = ismember(1:obj.beadData.numBeads,unique([obj.params.beadRangeToAnalyze.bead1 obj.params.beadRangeToAnalyze.bead2]));

            for bIdx =1:obj.beadData.numBeads
             obj.beadData.bead(bIdx).beadInRange = beadInRange(bIdx);    
            end
          
        end
        
        function ProcessDataBySegment(obj)
            % Find all unique bead lengths
            obj.segmentData.segmentLengthsInBeads = obj.beadData.allInteractions(:,2)-obj.beadData.allInteractions(:,1);
            % Create the encounter frequencies for each segment
            counter  = 1;
            for bIdx = 1:numel(obj.beadData.allInteractions(:,1))
                
                obj.segmentData.segment(counter).segmentRangeInBeads       = [obj.beadData.allInteractions(bIdx,1), obj.beadData.allInteractions(bIdx,2)];
                obj.segmentData.segment(counter).existInDB                 = true;
                obj.segmentData.segment(counter).segmentLengthInBeads      = obj.beadData.allInteractions(bIdx,2)- obj.beadData.allInteractions(bIdx,1);
                obj.segmentData.segment(counter).segmentPositionInTheChain = counter;
                
                % If there is a gap between the present segment and the next
                % one, fill it with demi segment
                if bIdx~=numel(obj.beadData.allInteractions(:,1)) && (obj.beadData.allInteractions(bIdx+1,1)-obj.beadData.allInteractions(bIdx,2))~=0
                    counter = counter+1;
                    obj.segmentData.segment(counter).segmentRangeInBeads       = [obj.beadData.allInteractions(bIdx,2), obj.beadData.allInteractions(bIdx+1,1)];
                    obj.segmentData.segment(counter).existInDB                 = false;
                    obj.segmentData.segment(counter).segmentLengthInBeads      = obj.beadData.allInteractions(bIdx+1,1)- obj.beadData.allInteractions(bIdx,2);
                    obj.segmentData.segment(counter).segmentPositionInTheChain = counter;
                end
                
                counter = counter+1;
            end
            
            obj.segmentData.histogram.freq   = hist([obj.segmentData.segment(:).segmentLengthInBeads],obj.beadData.numBeads);% length distribution
            obj.segmentData.histogram.bins   = obj.params.beadSizeInbp:obj.params.beadSizeInbp:obj.params.beadSizeInbp*obj.beadData.numBeads;
            obj.segmentData.histogram.xLabel = 'Estimated Segment Length';
            obj.segmentData.histogram.yLabel = 'Frequency';
            %            obj.PlotSegments
        end
        
        function [oneSide, twoSides,encounterOneSide,oneSideTr] = ProcessEncounters(obj,beadRangeToAnalyze,replicateName)
            % Process encounters             
            % Get encoutner matrix in the range specified by analysisRange
            allFreq   = obj.GetEncounterMatrix(beadRangeToAnalyze,replicateName);
            twoSides  = cell(numel(beadRangeToAnalyze.bead1),2);
            inds1     = beadRangeToAnalyze.bead1;
            inds2     = beadRangeToAnalyze.bead2;  
            encounterOneSide         = nan(obj.beadData.numBeads,obj.beadData.numBeads-1);
            
            % Compute  encounter probability 
            count = 1;
            for b1Idx = inds1
                for b2Idx = inds2
                    dist = abs(b1Idx-b2Idx);
                    if dist~=0
                    after  = nan;
                    before = nan;
                    if (b1Idx+dist)<=obj.beadData.numBeads-1
                        after  = allFreq(b1Idx,b1Idx+dist);
                    end
                    if (b1Idx-dist)>=1
                        before = allFreq(b1Idx,b1Idx-dist);
                    end
                    twoSides{count,1}(dist) = before;
                    twoSides{count,2}(dist) = after;                    
                    m    = obj.MeanIgnoreNaN( [before, after],2);
                    encounterOneSide(b1Idx,dist) = m;
                    end
                end
                % normalize two sides 
                s = obj.SumIgnoreNaN([twoSides{count,1},twoSides{count,2}],2);
                twoSides{count,1} = twoSides{count,1}/s;
                twoSides{count,2} = twoSides{count,2}/s;
                
                count = count+1;
            end
            
            % Define valid indices 
            oneSideTr = ~isnan(encounterOneSide);            
            oneSide   = cell(beadRangeToAnalyze.bead1(end)-beadRangeToAnalyze.bead1(1)+1,1);

       % Divide by the sum of each row to get probability 
            count = 1;
            for b1Idx= inds1
                encounterOneSide(b1Idx,:) = encounterOneSide(b1Idx,:)/obj.SumIgnoreNaN(encounterOneSide(b1Idx,:),2);
                encounterOneSide(b1Idx,~oneSideTr(b1Idx,:))= nan;
                oneSide{count,1} = encounterOneSide(b1Idx,:);
                count = count+1;
            end
        end
        
        function [encounterMatrix]   = GetEncounterMatrix(obj,analysisRange, replicateName)
            % Get encounter matrix of replicateNeame with data
            % corresponding to the analysis range. NaN is places in other
            % positions of the matrix.
            beadInds      = false(obj.beadData.numBeads);
            for b1Idx = analysisRange.bead1;
                for b2Idx = analysisRange.bead2
                    beadInds(b1Idx,b2Idx) = true;
                    beadInds(b2Idx,b1Idx) = true;
                end
            end
            encounterMatrix = obj.encounterMatrix.(lower(replicateName)).*double(beadInds);
            encounterMatrix(~beadInds) = nan;
        end
        
        function CreateEncounterFrequencyMatrices(obj)
            % For the two replicates and the average, create encounter
            % frequency matrices. These matrices represent the even
            % partition of the genome. Each bead's length is 3000bp.
            % This function calculate the encounter for all beads in the
            % database. 
            mumBeads                    = obj.beadData.numBeads;
            obj.encounterMatrix.rep1    = zeros(mumBeads);
            obj.encounterMatrix.rep2    = zeros(mumBeads);
            obj.encounterMatrix.average = zeros(mumBeads);
            
            for bIdx = 1:numel(obj.allData.beadStart1);
                bs1 = obj.allData.beadStart1(bIdx);%-obj.params.beadRangeToAnalyze(1)+1;
                be1 = obj.allData.beadEnd1(bIdx);%-obj.params.beadRangeToAnalyze(1)+1;
                bs2 = obj.allData.beadStart2(bIdx);%-obj.params.beadRangeToAnalyze(1)+1;
                be2 = obj.allData.beadEnd2(bIdx);%-obj.params.beadRangeToAnalyze(1)+1;
                for bsIdx = bs1:be1
                    for beIdx = bs2:be2
                        obj.encounterMatrix.rep1(bsIdx,beIdx)    = obj.encounterMatrix.rep1(bsIdx,beIdx)+obj.allData.rep1Count(bIdx);
                        obj.encounterMatrix.rep1(beIdx,bsIdx)    = obj.encounterMatrix.rep1(beIdx,bsIdx)+obj.allData.rep1Count(bIdx);
                        
                        obj.encounterMatrix.rep2(bsIdx,beIdx)    = obj.encounterMatrix.rep2(bsIdx,beIdx)+obj.allData.rep2Count(bIdx);
                        obj.encounterMatrix.rep2(beIdx,bsIdx)    = obj.encounterMatrix.rep2(beIdx,bsIdx)+obj.allData.rep2Count(bIdx);
                        
                        obj.encounterMatrix.average(bsIdx,beIdx) = obj.encounterMatrix.average(bsIdx,beIdx)+obj.allData.averageCount(bIdx);
                        obj.encounterMatrix.average(beIdx,bsIdx) = obj.encounterMatrix.average(beIdx,bsIdx)+obj.allData.averageCount(bIdx);
                    end
                end
            end
        end
        
        function FitData(obj)
            % Fit a function of the form y = a*x^-b to each monomer
            % encounter frequency data
            fNames = lower(obj.params.replicateName);% {'rep1','rep2','average'};
            model  = obj.params.fitModel; 
            for fIdx = 1:numel(lower(fNames));
                prevMissing = false;% indicator of whether the previous bead is missing.
                missingIdx  = [];% missing beads indices
                for bIdx = 1:numel(obj.params.beadRangeToAnalyze.bead1); %obj.beadData.numBeads
                    
                    freq = obj.beadData.encounterData.oneSide.(lower(fNames{fIdx}))(bIdx,1);% raw data
                    if ~all(freq{:}==0)% if it is present in the database
                        freq = {freq{:}./obj.SumIgnoreNaN(freq{:})};% normalize to create probability
                    end                    
                    y = freq{:}';
                    % Ignore missing data
                    includedPlaces = (y~=0 &~isnan(y));
                    y = y(includedPlaces);
                    x = find(includedPlaces);
                    
                    if ~all(y==0)
                        [fitObject, gof] = fit(x,y,model,obj.params.fOptions);
                        obj.beadData.bead(bIdx).fitResults.(fNames{fIdx})                 = obj.NewFitResultsStruct;
                        obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).bias            = (1/sum(x.^(-fitObject.slope)));% should divide by the actual distances 
                        obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).exp             = fitObject.slope;
                        obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).gof             = gof;
                        obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).beadDist        = x;
                        obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).encounterNumber = freq{:}(includedPlaces)';
                        obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).encounterProb   = y;%freq{:}';
                        obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).functionValues  = model(fitObject.slope,x);
                        obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).model           = model;
                        obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).existInDb       = true;                        
                        % save data in a struct for easy plotting
                        obj.results.fit.(fNames{fIdx}).exp(bIdx)                          = fitObject.slope;
                        obj.results.fit.(fNames{fIdx}).bias(bIdx)                         = obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).bias ;
                        
                        if prevMissing
                            if strcmpi(obj.params.fillGapsBy,'sameValuesAsBoundary')% options [sameValueAsBoundary/Interp]
                                % fill-in the gaps with fit values obtained at
                                % the boundary of the gap
                                if missingIdx(1)~=1
                                    bStart  = (obj.beadData.bead(missingIdx(1)-1).fitResults.(fNames{fIdx}).bias);
                                    bEnd    = (obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).bias);
                                    bValues = linspace(bStart,bEnd,numel(missingIdx));% linear interpolation of the bias values
                                    eStart  = obj.beadData.bead(missingIdx(1)-1).fitResults.(fNames{fIdx}).exp;
                                    eEnd    = obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).exp;
                                    eValues = linspace(eStart,eEnd,numel(missingIdx));% linear interpolation of the exp values
                                else
                                    % if the first data points are empty
                                    bValues = ones(1,numel(missingIdx))*(obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).bias);
                                    eValues = ones(1,numel(missingIdx))*(obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).exp);
                                end
                                
                                for mIdx = 1:numel(missingIdx)
                                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx})          = obj.NewFitResultsStruct;
                                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).bias     = bValues(mIdx);
                                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).exp      = eValues(mIdx);
                                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).gof      = 'missing bead- data interpolated from nearest neighbor';
                                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).beadDist = obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).beadDist;
                                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).encounterNumber = [];
                                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).encounterProb   = [];
                                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).functionValues  = [];
                                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).existInDb       = false;
                                    obj.results.fit.(fNames{fIdx}).exp(missingIdx(mIdx))  = obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).exp;
                                    obj.results.fit.(fNames{fIdx}).bias(missingIdx(mIdx)) = obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).bias;
                                end
                                missingIdx  = [];% reset the missing idx
                                prevMissing = false;% reset the missing flag
                            end
                        else
                            obj.results.fit.(fNames{fIdx}).exp(bIdx)  = fitObject.slope;
                            obj.results.fit.(fNames{fIdx}).bias(bIdx) = (1/sum((x).^(-fitObject.slope)));% fitObject.slope-1;
                        end
                    else
                        prevMissing = true;
                        missingIdx(end+1) = bIdx;% record the indices of the missing beads in the block
                    end
                end
                
                % Fill in the missing data of beads on the edge by the last
                % good observation                                
                for mIdx = 1:numel(missingIdx)
                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx})          = obj.NewFitResultsStruct;
                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).bias     = (obj.beadData.bead(missingIdx(1)-1).fitResults.(fNames{fIdx}).bias);
                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).exp      = obj.beadData.bead(missingIdx(1)-1).fitResults.(fNames{fIdx}).exp;
                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).gof      = 'missing bead- data interpolated from last observation';
                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).beadDist = obj.beadData.bead(missingIdx(1)-1).fitResults.(fNames{fIdx}).beadDist;
                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).encounterNumber = [];
                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).encounterProb   = [];
                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).functionValues  = [];
                    obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).existInDb       = false;
                    obj.results.fit.(fNames{fIdx}).exp(missingIdx(mIdx))  = obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).exp;
                    obj.results.fit.(fNames{fIdx}).bias(missingIdx(mIdx)) = obj.beadData.bead(missingIdx(mIdx)).fitResults.(fNames{fIdx}).bias;
                end
            end
        end
        
        function FitMeanModel(obj)
            % Calculate the mean of the mean data
            [distMin,distMax] = obj.GetDistanceRange;                 
            dists = distMin:distMax;            
            model                         = obj.params.fitModel;
            [meanSignal,~,~] = obj.GetMeanSignal(obj.params.beadRangeToAnalyze,'average');
             obj.results.fit.allData.meanEncounterProb = meanSignal;
            [fitObject, gof]              = fit(dists',meanSignal(dists)',model,obj.params.fOptions);
            obj.results.fit.allData.bias  = 1/(sum(dists.^(-fitObject.slope)));
            obj.results.fit.allData.exp   = fitObject.slope;
            obj.results.fit.allData.gof   = gof;
            obj.results.fit.allData.model = model;
            
            % revert weights 
            obj.params.fOptions.Weights      = [];
            
        end
        
        function PeakCalling(obj)
            % Find peaks in the one-sided encounter data
            % for TAD D and E seperately, and then analyze peaks between
            % TADs
            % The mean signal estimated by the model is used as the
            % expeected probability curve
            fNames     = lower(obj.params.replicateName);
            % for the case of 2 TADs , analyze each one seperately
            analysisRange(1).bead1 = 1:107;% analyze TAD D and the region between TADs 
            analysisRange(1).bead2 = 1:307;
            analysisRange(2).bead1 = 108:307;% analyze TAD E
            analysisRange(2).bead2 = 108:307;
%             analysisRange(3).bead1 = 1:107;
%             analysisRange(3).bead2 = 108:307;
            p = cell(numel(fNames), numel(analysisRange));
            for fIdx = 1:numel(fNames)
                obj.peaks.(lower(fNames{fIdx})) = [];
                for tIdx = 1:numel(analysisRange)
                    % Construct the one sided encounter matrix
                    [~, ~,encounterMat,~] = obj.ProcessEncounters(analysisRange(tIdx),fNames{fIdx});
                    % truncate the matrices 
                    startInd     = max(min(analysisRange(tIdx).bead2(1)-analysisRange(tIdx).bead1),1);
                    endInd       = max(analysisRange(tIdx).bead2(end)-analysisRange(tIdx).bead1);

                    encounterMat = encounterMat(analysisRange(tIdx).bead1,startInd:endInd);

                    obj.classes.peakFinder                = PeakCalling;
                    obj.classes.peakFinder.params.fitType = 'median';
                    obj.classes.peakFinder.FindPeaks(encounterMat);
                    % need to change peaks position according to analysis
                    % range 
                    
                    p{fIdx,tIdx} = obj.classes.peakFinder.peakList;       
%                      tempPeakList = p{fIdx,tIdx};
                    if ~isempty( p{fIdx,tIdx})
                     p{fIdx,tIdx}(:,1) =  p{fIdx,tIdx}(:,1)+analysisRange(tIdx).bead1(1)-1;
                     d          =   p{fIdx,tIdx}(:,2);% all distances
                     peaksHigh  = [p{fIdx,tIdx}(:,1), p{fIdx,tIdx}(:,1)+d];
                     indsHigh   = peaksHigh(:,2)>analysisRange(tIdx).bead2(end);
                     % Eliminate peaks over the max range
                     peaksHigh  = peaksHigh(~indsHigh,:); 
                     peaksLow   = [];                                 
                     peaksTotal = [peaksLow;peaksHigh];
                     
                     p{fIdx,tIdx} = peaksTotal;                                       
                    end
                    
%                     for pIdx = 1:size(p{fIdx,tIdx},1)
%                         tempPeakList(end+1,:) = p{fIdx,tIdx}(pIdx,:);
%                         
%                         tempPeakList(end,1) = tempPeakList(end,1)+analysisRange(tIdx).bead1(1)-1;
%                         % add peaks in both sides 
%                         d =  tempPeakList(end,2);
%                         if (tempPeakList(end,1)+d)<=(analysisRange(tIdx).bead2(end))
%                           tempPeakList(end,2)  =  tempPeakList(end,1)+d;                          
%                         end
%                         
%                         if ( tempPeakList(end,1)-d)>=analysisRange(tIdx).bead2(1) && tIdx~=3
%                              tempPeakList(end+1,1)   = tempPeakList(end,1);
%                              tempPeakList(end,2)     = tempPeakList(end,1)-d;
%                         end
%                         
%                     end
%                     p{fIdx,tIdx}= tempPeakList;
                end
                
                % keep only peaks in the analysis range 
                peaksList = cat(1,p{fIdx,:});
                r1    = obj.params.beadRangeToAnalyze.bead1;
                r2    = obj.params.beadRangeToAnalyze.bead2;
                pInds = (peaksList(:,1)>=r1(1) &  peaksList(:,2)<=r1(end)) & (peaksList(:,1)>=r2(1) &  peaksList(:,2)<=r2(end));
                obj.peaks.(fNames{fIdx}) = peaksList(pInds,:);
            end
        end
    
        function [meanSignal,meanSignalStd,tr] = GetMeanSignal(obj,analysisRange,fName)
            % tr are the valid index list when encounter data is sorte
            % according to distance [beads] 

            [distMin, distMax] = obj.GetDistanceRange;            
            dists = distMin:distMax;
            
            % Construct the valid index matrix 
            tr = false(obj.beadData.numBeads,obj.beadData.numBeads-1);
            for b1Idx = analysisRange.bead1;
                 if obj.beadData.bead(b1Idx).beadInRange
                  for b2Idx = analysisRange.bead2
                     if obj.beadData.bead(b2Idx).beadInRange       
                      tr(b1Idx, max(abs(b1Idx-b2Idx),1)) = true;
                      tr(b2Idx,max(abs(b1Idx-b2Idx),1)) = true;
                     end
                 end
                end
            end
            
            
            k = nan(numel(analysisRange.bead1),numel(dists));            
            for kIdx = 1:numel(analysisRange.bead1); % obj.beadData.numBeads-1
                if obj.beadData.bead(kIdx).beadInRange
                    inds         = find(tr(kIdx,:));
                    k(kIdx,inds) = obj.results.fit.(lower(fName)).bias(kIdx).*(inds).^(-obj.results.fit.(lower(fName)).exp(kIdx));
                    k(kIdx,~tr(kIdx,:))= nan;
                end
            end                     
                   
            meanSignal    = obj.MeanIgnoreNaN(k,1);
            meanSignalStd = obj.StdIgnoreNaN(k,1);

        end           
        
        function [oneSided,twoSided] = GetEncounterProbabilityMatrix(obj,analysisRange,fName)% nfinished
            numBeads = (obj.beadData.numBeads); 
            oneSided = zeros(numBeads,numBeads-1);             
            twoSided = zeros(numBeads,numBeads-1 +1 +numBeads-1);
            % The encounter signal probability
            for bIdx = 1:numel(analysisRange.bead1);%numBeads
                oneSided(bIdx,:) = obj.beadData.encounterData.oneSide.(lower(fName)){bIdx};
                t = obj.beadData.encounterData.twoSides.(fName)(bIdx,:);
                twoSided(bIdx,numBeads-numel(t{1}):numBeads-1) = fliplr(t{1});
                twoSided(bIdx,numBeads+1:numBeads+numel(t{2})) = t{2};
            end
        end
        
        function DisplayAllDataFit(obj,dispScale)
            % display mean encounter data fit
            if~exist('dispScale','var');
                dispScale = 'linear';
            end
            f = figure('FileName','meanDataFit');
            a = axes('Parent',f,...
                    'NextPlot','Add',...
                    'XScale',dispScale,...
                    'YScale',dispScale,...
                    'NextPlot','Add',...
                    'FontSize',40,...
                    'XLim',[1 obj.beadData.numBeads]);
            
            xlabel(sprintf('%s','Distance [beads]'),'FontSize',40);
            ylabel(sprintf('%s','Encounter Prob.'),'FontSize',40);
            title('Mean Data Fit','fontSize', 40)
            
            for bIdx = 1:numel(obj.params.beadRangeToAnalyze.bead1);
                if~isempty(obj.beadData.bead(bIdx).fitResults.average.encounterProb)
                    
                    lineC = rand(1,3);
                    line('XData',obj.beadData.bead(bIdx).fitResults.average.beadDist,...
                         'YData',obj.beadData.bead(bIdx).fitResults.average.encounterProb,...
                         'Color',lineC,...
                         'Marker','o',...
                         'MarkerSize',2,...
                         'MarkerEdgeColor',lineC,...
                         'MarkerFaceColor','k',...
                         'LineStyle','-',...
                         'LineWidth',3,...
                         'DisplayName',sprintf('%s%s','Bead',num2str(bIdx)),...
                         'handleVisibility','off',...
                         'Parent',a)
                else
                    %                         sprintf('%s%s%a', sprintf('%s%s',' fit for bead ',num2str((bIdx))),' is empty')
                end
            end
            
            % plot the fitted line
            A = obj.results.fit.allData.bias;
            B = obj.results.fit.allData.exp;
            [minDist,maxDist] = obj.GetDistanceRange;
            dists = minDist:maxDist;
            line('XData',dists,...
                'YData',(A*(dists).^(-B)),...
                'LineStyle','-',...
                'LineWidth',7,...
                'Color','r',...
                'DisplayName',['mean data fit. \beta =' num2str(B)] );
            legend(get(a,'Children'));
            
        end
        
        function DisplayEncounterMatrices(obj, windowSize)
            % Display the encounter matrices. the pixels of the matrices
            % indicate the number of times bead i met bead j.
            % the parameter windowsize, determine the window size of the
            % median filter. The median filter is utilized to follow the
            % protocol used by Nora et al 2012.
            % windowSize must be a positive integer.
            if ~exist('windowSize','var')
                windowSize = 1;
            end
            cmap   = hot(256);% define color map
            fNames = (obj.params.replicateName);
            titles = {'Replicate 1', 'Replicate 2', 'Average'};
            [minDist,maxDist] = obj.GetDistanceRange;
            dists = minDist:maxDist;
            for fIdx = 1:numel(fNames)
                figure('Name',sprintf('%s%s','Encounter Matrix ',fNames{fIdx}),'FileName',sprintf('%s%s', 'encounterMatrix',fNames{fIdx}));
                em    = obj.GetEncounterMatrix(obj.params.beadRangeToAnalyze, lower(fNames{fIdx}));
                em = em(obj.params.beadRangeToAnalyze.bead1,dists);
                if windowSize>1
                em    = medfilt2(em,[windowSize,windowSize]);
                
                maxEm = max(em(:));
                em    = em./maxEm;% normalize
                end
                % construct a color scheme from white to black through red such
                % that black is
                imagesc(em);
                set(gca,'FontSize',40)
                title(gca,titles{fIdx},'FontSize',40);
                colormap((cmap))
                axis ij
            end

        end
        
        function DisplayEncounterProbabilityByBead(obj,beads,dispScale,sides,distance)
            % Plot the encounter probability curves for the beads indicated
            % in the numerical array beads
            % with a display scale dispScale ='linear'/'log'
            % display sides= 'oneSided'/'twoSides'
            % at a distance indicated by the numerical array distance
            
            if ~exist('beads','var')
                beads = obj.params.beadRangeToAnalyze.bead1;% 1:numel(obj.beadData.bead);
            end
            if ~exist('dispScale','var')
                dispScale = 'linear';% axes display scale [linear/log]
            end
            
            if ~exist('sides','var')
                sides = 'oneSide';
            end
            

            fNames = obj.params.replicateName;
            for fIdx = 1:numel(fNames)
                % create main figure
                f = figure('FileName',['encounterProbabilityByDistance',fNames{fIdx}],...
                    'Name',['encounterProbabilityByDistance',fNames{fIdx}]);
                [minDist,maxDist] = obj.GetDistanceRange;
                if strcmpi(sides,'oneSide')
                    xLim = [minDist,maxDist]; 
                elseif strcmpi(sides,'twoSides')
                    xLim = [-maxDist,maxDist];
                end
                % create main axes
                a = axes('Parent',f,...
                        'NextPlot','Add',...
                        'XScale',dispScale,...
                        'YScale',dispScale,...
                        'FontSize',40,...
                        'XLim',xLim);
                
                xlabel(sprintf('%s','Distance [beads]'),'FontSize',40);
                ylabel(sprintf('%s','Prob. Encounter'),'FontSize',40);
                
                title(fNames{fIdx},'fontSize',40)
                for bIdx = beads
                    if~isempty(obj.beadData.bead(bIdx).fitResults.(lower(fNames{fIdx})).encounterProb)
                        fResults = obj.beadData.bead(bIdx).fitResults.(lower(fNames{fIdx}));
                        % indicate whether the bead is missing in the
                        % legend
                        if fResults.existInDb
                            dispName = sprintf('%s%s','Bead',num2str(bIdx));
                        else
                            dispName = sprintf('%s%s%s','Bead',num2str(bIdx),' (missing)');
                        end
                        
                        % plot one side encounter prob
                        if strcmpi(sides,'oneSide')
                          if ~exist('distance','var')
                              distance = fResults.beadDist; % at what distance to display the encounter probability 
                          end
                            inds = distance(distance<=numel(fResults.encounterProb));
                            line('XData',inds,...
                                'YData',fResults.encounterProb(inds),...
                                'Color',rand(1,3),...
                                'Marker','.',...
                                'MarkerSize',7,...
                                'LineStyle','-',...
                                'LineWidth',3,...
                                'DisplayName',dispName,...
                                'Parent',a);
                            
                           
                        elseif strcmpi(sides,'twoSides')
%                               eData = [fliplr(obj.beadData.encounterData.twoSides.(lower(fNames{fIdx})){bIdx,1}),...
%                                      obj.beadData.encounterData.twoSides.(lower(fNames{fIdx})){bIdx,2}];                                 
                            if ~exist('distance','var')
                              distance = (1:numel(obj.beadData.encounterData.twoSides.(lower(fNames{fIdx})){bIdx,2})); % at what distance to display the encounter probability 
                              indsLeft  = distance;
                              indsRight = distance;
                            else
                               indsLeft  = distance;
                               indsRight = distance;
                            end
                           
                            eData     = [fliplr(obj.beadData.encounterData.twoSides.(lower(fNames{fIdx})){bIdx,1}(indsLeft)),...
                                         obj.beadData.encounterData.twoSides.(lower(fNames{fIdx})){bIdx,2}(indsRight)];
                                                
                            line('XData',[-fliplr(indsLeft) indsRight],...
                                'YData',eData,...
                                'Color',rand(1,3),...
                                'Marker','.',...
                                'MarkerSize',7,...
                                'LineStyle','-',...
                                'LineWidth',3,...
                                'DisplayName',dispName,...
                                'Parent',a);            
                             clear distance
                        end
                    else
                        % do nothing
                    end
                    
                end
            end
        end
        
        function DisplayEncounterProbabilityByPosition(obj,beadRange,pos,fName)
            %%% Display the encounter probability at a specific position
            % beadRang is a n integer list (sorted) for the beads for which
            % to display the encounter probability 
        [oneSided,~] = obj.GetEncounterProbabilityMatrix(obj.params.beadRangeToAnalyze,fName);
        f = figure('Name','EncounterProbByDistance',...
                   'Units','norm');
        a = axes('Parent',f,...
                 'Units','norm',...
                 'FontSize',40,...
                 'LineWidth',4,...
                 'NextPlot','Add');
        title(a,['Encounter Prob by bead, distance ' num2str(pos)],'FontSize',40);
         ep = oneSided(beadRange,pos);
         line('XData',beadRange,...
              'YData',ep,...
              'LineStyle','none',...
              'Marker','o',...
              'MarkerSize',9,...
              'MarkerFaceColor','b',...
              'Parent',a,...
              'DisplayName','Encounter Prob');
         xlabel(a,'Bead number','FontSize',40);
         ylabel(a,['P(' num2str(pos) ')'],'FontSize',40);
         
         % Insert the mean (ignore nan and zeros)
         m = obj.MeanIgnoreNaN(ep(ep>0)); % do not include zero observations
         plot(a,[1 beadRange(end)],[m m],'-.r',...
             'LineWidth',5,'DisplayName',['Mean=' num2str(m)]);
         legend(get(a,'Children'))
        end
        
        function DisplayEncounterProbabilityHistogramByPosition(obj,beadRange,pos,fName)
          [oneSided,~] = obj.GetEncounterProbabilityMatrix(obj.params.beadRangeToAnalyze,fName);
          ep =  oneSided(beadRange,pos);
          ep = ep(ep>0);
          ep = ep(~isnan(ep));
          [h,b] = hist(ep,floor(numel(beadRange)/2));
          f = figure('Name','probability histogram');
          a = axes('Parent',f,...
                   'Units','norm',...
                   'FontSize',40,...
                   'LineWidth',4);
          bar(a,b,h)
          xlabel(a,['P(', num2str(pos),')'],'FontSize',40)
          ylabel(a,'Freq','FontSize',40);
          
               
        end
        
        function DisplayFittedParameters(obj)
            % Display the fitted parameters of the encounter probability, using the model defined in obj.modelfit,
            % for the two replicas and their average
            fNames = obj.params.replicateName;% {'Rep1','Rep2','Average'};
            
            for fIdx = 1:numel(fNames)
                % plot fitted beta values
                expFigName = sprintf('%s%s','\beta for ', fNames{fIdx});
                
                % create main figure
                fe = figure('Name',expFigName,...
                            'FileName',['FittedExpValues',fNames{fIdx}]);
                % create main axes
                ae = axes('Parent',fe,...
                          'FontSize',40,...
                          'XLim',[obj.params.beadRangeToAnalyze.bead1(1) obj.params.beadRangeToAnalyze.bead1(end)],...
                          'NextPlot','Add');
                
                title(ae,expFigName,'FontSize',40);
                xlabel(ae,'Bead number','FontSize',40);
                ylabel(ae,'Fitted \beta','FontSize',40);
                
                % plot fitted bias values
                biasFigName = sprintf('%s%s','Bias for ', fNames{fIdx});
                
                fb = figure('Name',biasFigName,...
                            'FileName',['FittedExpValues',fNames{fIdx}]);
                ab = axes('Parent',fb,...
                         'FontSize',40,...
                         'XLim',[obj.params.beadRangeToAnalyze.bead1(1) obj.params.beadRangeToAnalyze.bead1(end)],...
                         'NextPlot','Add');

                title(ab,biasFigName,'FontSize',40);
                xlabel(ab,'Bead number','FontSize',40);
                ylabel(ab,'Fitted bias','FontSize',40);
                
                
                for bIdx = 1:numel(obj.params.beadRangeToAnalyze.bead1);
                    % indicate whether the bead is missing in the legend
                    if obj.beadData.bead(bIdx).fitResults.(lower(fNames{fIdx})).existInDb
                        dispName = ['bead ', num2str(obj.params.beadRangeToAnalyze.bead1(bIdx))];
                    else
                        dispName = ['bead ', num2str(obj.params.beadRangeToAnalyze.bead1(bIdx)),' (missing)'];
                    end
                    % plot the values of the fitted exponent
                    line('Parent',ae,...
                        'XData',obj.params.beadRangeToAnalyze.bead1(bIdx),...
                        'YData',obj.beadData.bead(bIdx).fitResults.(lower(fNames{fIdx})).exp,...
                        'LineStyle','none',...
                        'Marker','o',...
                        'MarkerFaceColor','b',...
                        'DisplayName',dispName)
                    
                    % plot the values of the fitted bias
                    line('Parent',ab,...
                        'XData',obj.params.beadRangeToAnalyze.bead1(bIdx),...
                        'YData',obj.beadData.bead(bIdx).fitResults.(lower(fNames{fIdx})).bias,...
                        'Marker','o',...
                        'MarkerFaceColor','b',...
                        'LineStyle','none',...
                        'DisplayName',dispName);
                end
                
            end
        end
        
        function DisplayFitByBead(obj,beadNumbers,dispScale)
            warning off
            fNames = lower(obj.params.replicateName);% {'rep1','rep2','average'};
            if ~exist('beadNumbers','var')||strcmpi(beadNumbers,'all')
                beadNumbers = obj.params.beadRangeToAnalyze.bead1;
            end
            
            if ~exist('dispScale','var')
                dispScale = 'linear';
            end
            [minDist, maxDist] = obj.GetDistanceRange;
            
            for fIdx = 1:numel(fNames)
                f = figure('Name',['FitByBead',fNames{fIdx}]);
                a = axes('Parent',f,...
                        'NextPlot','Add',...
                        'XLim',[minDist maxDist],...
                        'XScale',dispScale,...
                        'YScale',dispScale,...
                        'NextPlot','Add',...
                        'FontSize',40);
                
                xlabel(sprintf('%s','Distance [beads]'),'FontSize',40);
                ylabel(sprintf('%s','Num. encounters'),'FontSize',40);
                title(fNames{fIdx},'FontSize',40)
                % line color
                lineC = [linspace(0,1,numel(beadNumbers))',0.5*ones(numel(beadNumbers),1),0.5*ones(numel(beadNumbers),1)];
                
                for bIdx = 1:numel(beadNumbers)
                    if~isempty(obj.beadData.bead(beadNumbers(bIdx)).fitResults)
                        %                         freq = obj.beadData.bead(beadNumbers(bIdx)).fitResults.(fNames{fIdx}).encounterNumber;
                        %                         freq(freq==0)=NaN;
                        % indicate in the legend if the bead is missing
                        if obj.beadData.bead(beadNumbers(bIdx)).fitResults.(fNames{fIdx}).existInDb
                            dispName = sprintf('%s%s','Bead',num2str(beadNumbers(bIdx)));
                        else
                            dispName = sprintf('%s%s','Bead',num2str(beadNumbers(bIdx)),' (missing)');
                        end
                        line('XData',obj.beadData.bead(beadNumbers(bIdx)).fitResults.(fNames{fIdx}).beadDist,...
                            'YData',obj.beadData.bead(beadNumbers(bIdx)).fitResults.(fNames{fIdx}).encounterProb,...
                            'Color',lineC(bIdx,:),...
                            'Marker','none',...
                            'MarkerSize',2,...
                            'MarkerEdgeColor','c',...
                            'MarkerFaceColor','k',...
                            'LineStyle','-',...
                            'LineWidth',2,...
                            'DisplayName',dispName,...
                            'Parent',a)
                        
                        % plot the fitted line
                        A = obj.beadData.bead(beadNumbers(bIdx)).fitResults.(fNames{fIdx}).bias;
                        B = obj.beadData.bead(beadNumbers(bIdx)).fitResults.(fNames{fIdx}).exp;
                        line('XData',obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).beadDist,...
                            'YData',(A*(obj.beadData.bead(bIdx).fitResults.(fNames{fIdx}).beadDist).^(-B)),...
                            'LineStyle','-',...
                            'LineWidth',4,...
                            'Color','r',...
                            'DisplayName',['Bead ',num2str(beadNumbers(bIdx)),' \beta: ',num2str(B), ' bias: ',num2str(A)]);
                        
                    else
                        sprintf('%s%s%a', sprintf('%s%s','Bead',num2str(beadNumbers(bIdx))),'is empty');
                    end
                end
            end
        end
        
        function DisplayPeaks(obj)
            % Display binary matrix of peaks             
            fNames = lower(obj.params.replicateName);
            for fIdx = 1:numel(fNames)
                
                f = figure('Name',sprintf('%s%s', 'Peak List ' , fNames{fIdx}),...
                           'FileName',sprintf('%s%s', 'Peak List' , fNames{fIdx}),...
                           'Units','norm');
                a = axes('Parent',f,...
                         'Units','norm',...
                         'NextPlot','Add',...
                         'FontSize',40,...
                         'LineWidth',4);
                 
                 
                 e = obj.encounterMatrix.(fNames{fIdx});
                 e = e./max(e(:));                                  
                 s = surface(e,'EdgeColor','none',...
                               'FaceLighting','phong',...
                               'FaceColor','interp');
                 % add light 
                 l = light('Parent',a,...
                           'Position',[-0.2 0.1 1],...
                           'HandleVisibility','off');
                 camlight left
                  plot3(a,obj.peaks.(fNames{fIdx})(:,1),obj.peaks.(fNames{fIdx})(:,2),...
                        e(obj.peaks.(fNames{fIdx})(:,1),...
                          obj.peaks.(fNames{fIdx})(:,2)),'or','MarkerSize',8)

                     title(a,fNames{fIdx});
                     xlabel(a,'Bead');
                     ylabel(a,'Bead');
                     colormap summer
                     
            end
        end
        
        function PlotSegments(obj)
            % plot a cylinder representing the polymer, with red
            % representing missing segments in the data
            maxBead     = obj.beadData.numBeads;
            allSegments = obj.beadData.allInteractions;
            
            segmentLength = 1;
            pLength       = maxBead;
            
            % plot the missing pieces
            f = figure('MenuBar','none');
            a = axes('Parent',f,...
                'NextPlot','Add');
            % create a cylinder lines
            t  = linspace(0,2*pi,50);
            r  = 10;
            x  = repmat(r*cos(t),pLength,1);
            y  = repmat(r*sin(t),pLength,1);
            z  = repmat((1:pLength)',1,numel(t));
            c  = 255*ones(size(x));
            m  = mesh(x,y,z,'Parent',a,...
                'CData',c,...
                'FaceColor','b',...
                'CDataMapping','scaled',...
                'FaceLighting','gouraud');
            
            for bIdx = 1:numel(allSegments(:,1))
                % plot places where we know segments exist for the
                % experiments
                c   = get(m,'CData');
                c((allSegments(bIdx,1)-1)*segmentLength+1:(allSegments(bIdx,2))*segmentLength,:)=0;
                set(m,'CData',c);
            end
            
            daspect([1 1 1]);
            cameratoolbar
        end
        
        function obj = ReloadObj(obj,structIn)
            % Reload a constructed obj with a structure
            prop= properties(structIn);
            for fIdx = 1:numel(prop)
                obj.(prop{fIdx}) = structIn.(prop{fIdx});
            end
        end
        
        function objOut = SaveObj(obj)
            % Save the object to a structure for saving as mat file
            prop= properties(obj);
            for fIdx = 1:numel(prop)
                objOut.(prop{fIdx}) = obj.(prop{fIdx});
            end
            uisave('objOut')
        end
        
        function FitCurveToExpData(obj,smoothingFac)
            % fit a curve to the exponents of the fitted encounter
            % frequency data. use a smoothing spline with a smoothingFac>0
            if ~exist('smoothingFac','var')
                smoothingFac = 1;
            end
            fNames = obj.params.replicateName;% {'Rep1','Rep2','Average'};
            % fit using a smoothing spline
            for fIdx = 1:numel(fNames)
                [obj.results.fit.(lower(fNames{fIdx})).expSpline,~] = spaps(1:numel(obj.results.fit.(lower(fNames{fIdx})).exp),...
                    obj.results.fit.(lower(fNames{fIdx})).exp',smoothingFac);
                
                f = figure('Name',['fittedExpValuesWithSpline',fNames{fIdx}],...
                    'FileName',['fittedExpValuesWithSpline',fNames{fIdx}]);
                a = axes('Parent',f,...
                    'NextPlot','Add',...
                    'FontSize',40,...
                    'XLim',[1 obj.beadData.numBeads],...
                    'YLim',[0 max(obj.results.fit.(lower(fNames{fIdx})).exp)+0.2]);
                
                line('XData',1:numel(obj.beadData.bead),...
                    'YData',obj.results.fit.(lower(fNames{fIdx})).exp,...
                    'Marker','o',...
                    'MarkerEdgeColor','k',...
                    'MarkerSize',8,...
                    'MarkerFaceColor','b',...
                    'lineStyle','none',...
                    'DisplayName','\beta',...
                    'Parent',a);
                
                title(sprintf('%s%s','Fitted \beta for ', fNames{fIdx}),'FontSize',40);
                xlabel('Bead number','FontSize',40);
                ylabel('Fitted \beta values','FontSize',40);
                
                % plot the smoothing spline
                fn = fnplt(obj.results.fit.(lower(fNames{fIdx})).expSpline);
                line('XData',fn(1,:),...
                    'YData',fn(2,:),...
                    'DisplayName','Spline',...
                    'Parent',a,...
                    'LineWidth',3,...
                    'Color','k');
                
                legend(get(a,'Children'))
            end
        end
        
        function CheckEncounterDataSymmetry(obj)
            % Check the two-sided encounter data for each bead.
            % The encounter data of 'left' and 'right' are compared by
            % taking the difference of the minimal number of beads of both
            % sides.
            fNames = obj.params.replicateName;% {'Rep1','Rep2','Average'};
            d      = zeros(obj.beadData.numBeads,numel(fNames));
            meanDiff = zeros(1,numel(fNames));
            for fIdx = 1:numel(fNames)
                for bIdx = 1:numel(obj.beadData.encounterData.twoSides.(lower(fNames{fIdx}))(:,1))
                    % find the end with the least number of beads
                    [~,m] = min([numel(obj.beadData.encounterData.twoSides.(lower(fNames{fIdx})){bIdx,1}),...
                        numel(obj.beadData.encounterData.twoSides.(lower(fNames{fIdx})){bIdx,2})]);
                    numPoints = numel(obj.beadData.encounterData.twoSides.(lower(fNames{fIdx})){bIdx,m});
                    
                    if ~isempty(obj.beadData.encounterData.twoSides.(lower(fNames{fIdx})){bIdx,1})
                        
                        pLeft = obj.beadData.encounterData.twoSides.(lower(fNames{fIdx})){bIdx,1}(1:numPoints);
                        pLeft = pLeft./obj.SumIgnoreNaN(pLeft);
                        
                        pRight = obj.beadData.encounterData.twoSides.(lower(fNames{fIdx})){bIdx,2}(1:numPoints);
                        pRight = pRight./obj.SumIgnoreNaN(pRight);% get prob by normalizeing
                        
                        d(bIdx,fIdx) = mean(pLeft-pRight);
                    else
                        d(bIdx,fIdx) = NaN;
                    end
                end
                
                meanDiff(fIdx) = mean(d(~isnan(d(:,fIdx)),fIdx));
                
                f(fIdx) = figure;
                ax(fIdx) = axes('Parent',f(fIdx),'NextPlot','Add');
                
                % plot the mean prob difference of the two sides for each bead
                line('XData',1:obj.beadData.numBeads,...
                    'YData',d(:,fIdx),'Linewidth',6,...
                    'Parent',ax(fIdx),...
                    'DisplayName','two-sided encounter Prob diff.',...
                    'Color','b')
                % plot the mean prob difference of the two sides for all beads
                line('XData',[1,numel(d(:,fIdx))],...
                    'YData',meanDiff(fIdx)*ones(1,2),...
                    'Color','r',...
                    'LineWidth',6,...
                    'Parent',ax(fIdx),...
                    'DisplayName','mean prob diff')
                
                xlabel('Bead number','FontSize',40);
                ylabel('mean encounter prob. difference','FontSize',40);
                title(fNames{fIdx},'Fontsize',40)
                set(gca,'Fontsize',40)
                legend(get(ax(fIdx),'Children'));
            end
            
        end
    end
    
    methods (Static)
        
        function objOut = LoadObj(objIn)
            % Constuct the objec with a structure objIn
            if nargin<1
                objIn = uiload('*.mat');
            end
            objOut = AnalyzeEncounterFrequencies;
            objOut = ReloadObj(objOut,objIn);
        end
        
        function fitRStruct = NewFitResultsStruct()
            % create a new fitting structure to hold results
            fitRStruct = struct('bias',[],...
                'exp',[],...
                'gof',[],...
                'beadDist',[],...
                'encounterNumber',[],...
                'encounterProb', [],...
                'functionValues',[],...
                'existInDb',[],...              
                'model',[]);
        end
        
        function valsOut = model(x,encounterData)%unused
            % this ia the SSD of the fitted parameters
            % encounterData is a two column vector. the first column is the
            % encounter frequencies, the second is the index of included
            % beads.
            % encounterData(:,1) - log(encounter frequencies)
            % encounterData(:,2) - log(included bead indices)
            % x(1) - is the bias
            % x(2) - is the slope
            logBeadDist      = encounterData(:,2);
            logEncounterFreq = encounterData(:,1);
            valsOut = (sum((x(1)-x(2)*logBeadDist-logEncounterFreq).^2));
        end
        
        function [valInEq, valEq] = modelConstraint(x,includedBeads)
            % constrain of the fitted model
            % we wish the sum of elements to be one
            logBeadDist = includedBeads;
            valEq       = sum(exp(x(1)-x(2)*logBeadDist))-1;% *numel(includedBeads)-x(2)*sum((includedBeads))-1;
            valInEq     = valEq;
        end
        
        function sigOut = LoessSmooth(x,y,span)%unused
            % Smooth the signal sigIn with a loess filter
            sigOut = smooth(x,y,span,'loess');
        end
        
        function m = MeanIgnoreNaN(vec,dim)
            if ~exist('dim','var')
                dim = 1;
            end
            
            m = AnalyzeEncounterFrequencies.SumIgnoreNaN(vec,dim);
            m = m./sum(~isnan(vec),dim);
                        
        end
        
        function s = SumIgnoreNaN(vec,dim)
            % sum of matrix or vaector ignoring NaN values (replacing them
            % with zeros)
            if ~exist('dim','var')
                if size(vec,1)==1                
                    dim = 2;
                else
                    dim = 1;
                end
            end
            
            n = isnan(vec);
            vec(n) = 0;
            s = sum(vec,dim);            
        end
        
        function st = StdIgnoreNaN(vec,dim)
            ind = ~isnan(vec);
            if ~exist('dim','var')
                dim = 1;
            end                        
                
            for dIdx = 1:size(vec,2) 
                if dim==1
                    sig   = vec(ind(:,dIdx),dIdx);
                    st(1,dIdx) = std(sig);
                else
                    sig = vec(dIdx,ind(dIdx,:));
                    st(dIdx,1) = std(sig);
                end
              
            end
        end
            
    end
    
    % private methods
    methods (Access=private)
        
        function CheckAnalysisRangeMonotonicity(obj, analysisRange)
            % Produce an error if the analysis range is non-monotonically
            % increasing 
            d1 = diff(analysisRange.bead1);
            d2 = diff(analysisRange.bead2);
            if any(d1<0) || any(d2<0) 
               error('Analysis range must be a monotonic vector of integers') 
            end
        end
        
        function [distMin, distMax] = GetDistanceRange(obj)
            % Get the minimum and maximal range of distances between
            % beads in the analysis range 
            a1 = obj.params.beadRangeToAnalyze.bead1(1);
            b1 = obj.params.beadRangeToAnalyze.bead1(end);
            a2 = obj.params.beadRangeToAnalyze.bead2(1);
            b2 = obj.params.beadRangeToAnalyze.bead2(end);
            distMax = max([abs(b2-a1), abs(b1-a2)]);          
            distMin = min([[abs(a1-b2),abs(a2-b1)],max(min([b2-a1,a2-b1]),1)]);
        end
    end
    
end