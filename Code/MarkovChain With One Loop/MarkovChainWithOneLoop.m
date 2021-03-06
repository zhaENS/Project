classdef MarkovChainWithOneLoop<handle
    % This function represent the solution of the markov chain representing a
    % chain with a single loop. one end is loose and can connect to any other
    % bead to create a loop (excluding itself).
    properties
        params
        transitionMatrix
        solution
        fitResults
    end
    
    properties (Access=private)
        sol      = @(t,p0,v,e)(v)*expm(e*t)/(v)*p0;
        modelFit = fittype('(1/sum([1:numel(x)].^(-b))).*x.^(-b)');
    end
    
    methods
        
        function  obj = MarkovChainWithOneLoop
            obj.Initialize;
        end
        
        function DefineModelParams(obj)
            % Define the parameters for the chain
            obj.params.numBeads        = 64;
            obj.params.encounterDist   = 0.1; % the maximal distance two bead are considered to be in contact
            obj.params.connectorLength = 0.1;
            obj.params.b               = 1;  % dist STD between beads 
            % calculate the attachment rate to any other bead (will be changed according to distance) when free
            
%             % the PDF of the bead distance
%             f = @(r,n,b)((3/(2*pi*n*b^2))^(3/2)).*(exp(-((3*r.^2)./(2*n*b^2))));
%             t = 0:0.001:obj.params.encounterDist;
%             t = linspace(-2*obj.params.numBeads*obj.params.connectorLength,2*obj.params.numBeads*obj.params.connectorLength,1001);
            
            %             for  aIdx = 1:obj.params.numBeads-1
            %
            % %                 obj.params.alpha(aIdx) = ((3/(2*pi*obj.params.connectorLength^2 * aIdx))^(3/2))*trapz(t,exp(-3*(t.^2)./(2*aIdx*obj.params.connectorLength^2)));
            %                   obj.params.alpha(aIdx) = f(t,aIdx,obj.params.connectorLength);
            %                   F = f(t,1,0.1);
            %                 tF = cumsum(F);
            %                 cF = tF./max(tF(:));
            %             end
            
            
            %              obj.params.alpha           = normcdf(obj.params.encounterDist,0,sqrt((1:obj.params.numBeads-1))*obj.params.connectorLength)-...
            %                                           normcdf(sqrt(3)*obj.params.encounterDist,0,0);
            %             obj.params.alpha           = 0.1*ones(1,obj.params.numBeads-1);
            k                          = 2:obj.params.numBeads;
            normConst = sum(((1./abs(1-k)).^1.5).*exp(-3*(obj.params.encounterDist.^2)./(2* obj.params.b^2 *abs(1-k))));
            obj.params.alpha           = (1/normConst)*((1./abs(1-k)).^1.5).*exp(-3*(obj.params.encounterDist.^2)./(2* obj.params.b^2 *abs(1-k)));
%             obj.params.alpha           = (1/sum((1:obj.params.numBeads-1).^(-1.5))).*((1:obj.params.numBeads-1).^(-1.5));
            obj.params.beta            = 1;%-obj.params.alpha;%ones(1,obj.params.numBeads-1)*(1/3);% the detachment rate
            obj.params.mu              = ones(1,obj.params.numBeads-1)*(1/3); % rate of forming a smaller loop when attached (with the neighbor bead)
            obj.params.lambda          = ones(1,obj.params.numBeads-1)*(1/3); % rate of forming a bigger loop when attached (with the neighbor bead)
            obj.params.timePoints      = 0:0.01:(25-0.01);
            obj.params.connectedBeads  = [];
            
            % deine beads which are hot spots and their transition rates
            obj.params.hotSpots       = [];
            obj.params.hotSpotMu      = [];%obj.params.mu(1)/4;
            obj.params.hotSpotLambda  = [];%obj.params.lambda(1)/4;
            obj.params.hotSpotBeta    = [];%obj.params.beta(1)/4;
            % TODO: introduce attachment distance into the alpha parameter
            obj.params.initialDistribution = rand(1,obj.params.numBeads);% (1/obj.params.numBeads)*ones(obj.params.numBeads,1); % initial probability of states
            obj.params.initialDistribution =obj.params.initialDistribution'/sum(obj.params.initialDistribution );
            
            for hIdx = 1:numel(obj.params.hotSpots)
                obj.params.beta(obj.params.hotSpots(hIdx))     = obj.params.hotSpotBeta;
                obj.params.lambda(obj.params.hotSpots(hIdx)+1) = obj.params.hotSpotLambda;
                obj.params.mu(obj.params.hotSpots(hIdx))       = obj.params.hotSpotMu;
            end
            
        end
        
        function Initialize(obj)
            % Create the transition probability matrix
            
            obj.DefineModelParams
            
            obj.CreateTransitionMatrix
            
            obj.Solve
            
            obj.PlotProbabilityEvolution
            
            obj.FitEncounterProbability
        end
        
        function CreateTransitionMatrix(obj)
            %             % create the transition matrix
            %             M = zeros(obj.params.numBeads);
            %             M = M+diag(obj.params.lambda,-1);
            %             M = M+diag(obj.params.mu,1);
            %             M(1,2:end) = obj.params.beta;
            %             M(2:end,1) = obj.params.alpha;
            %             M = M-diag(sum(M));
            %             obj.transitionMatrix = M;
            
            M          = zeros(obj.params.numBeads);
            M          = M - diag([1, (obj.params.beta+obj.params.lambda+obj.params.mu)]);% set main diagonal
            M(1,1)     = -1;
            M(2,2)     = -(obj.params.lambda(1)+obj.params.beta(1));
            M(end,end) = -(obj.params.mu(1)+obj.params.beta(1));
            M          = M+diag(obj.params.mu,1); % set super diagonal
            M          = M+diag(obj.params.lambda,-1); % set sub-diagonal
%             M(2:end,1) = obj.params.alpha;%(1./sum((1:obj.params.numBeads-1).^(-1.5)))*(1:obj.params.numBeads-1)'.^(-1.5);% set first row alpha(k)
%             M(1,2:end) = obj.params.beta;% set first row;
            
            obj.transitionMatrix = M;
        end
        
        function Solve(obj)
            % solve the differential equations of the chain 
            [v,D] = eig(obj.transitionMatrix);
            obj.solution = zeros(obj.params.numBeads,numel(obj.params.timePoints));
            for tIdx = 1:numel(obj.params.timePoints)
                obj.solution(:,tIdx)= obj.sol(obj.params.timePoints(tIdx),obj.params.initialDistribution,v,D);
            end
        end
        
        function PlotProbabilityEvolution(obj)
            f = figure('Name','MarkovChainWithOneLoopProbabilityAtStSt');
            a = axes('Parent',f,...
                     'NextPlot','Add',...
                     'FontSize',40);            
            for bIdx = 1:obj.params.numBeads
                dName = ['bead',num2str(bIdx)];
                if bIdx ==1
                    dName = 'no loop';
                end
                line('XData',obj.params.timePoints,...
                     'YData',obj.solution(bIdx,:),...
                     'Color',rand(1,3),...
                     'DisplayName',dName,...
                     'LineWidth',5,...
                     'Parent',a);
            end
            
            xlabel(a,'Time','FontSize',40);
            ylabel(a,'Prob.','FontSize',40);            
            legend(get(a,'Children'));
           
        end
        
        function FitEncounterProbability(obj)
            % Fit a curve to the steady state probability of loops
            [obj.fitResults.fitParams,obj.fitResults.gof] = fit((1:obj.params.numBeads-1)',obj.solution(2:end,end),obj.modelFit,...
                'StartPoint',1.5,...
                'robust','LAR',...
                'TolFun',1e-16,...
                'TolX',1e-16);
            f = figure('Name','Loop Probability At Steady State');
            a = axes('Parent',f,'NextPlot','Add','FontSize',40);
            
            % draw the stst probabilites 
            line('XData',1:obj.params.numBeads-1,...
                 'YData',obj.solution(2:end,end),...
                 'Color','b',...
                 'LineWidth',5,...
                 'Parent',a,...
                 'DisplayName','Solution at StSt');
             
             % draw the fitted curve 
            line('XData',1:obj.params.numBeads-1,...
                 'YData',(1/sum((1:obj.params.numBeads).^(-obj.fitResults.fitParams.b)))*(1:obj.params.numBeads-1).^(-obj.fitResults.fitParams.b),...
                 'Color','r',...
                 'LineWidth',5,...
                 'Parent',a,...
                 'DisplayName',['fitted \beta = ' num2str(obj.fitResults.fitParams.b)]);
             
            % draw the expected encounter probability in a free chain
             line('XData',1:obj.params.numBeads-1,...
                  'YData',(1/sum((1:obj.params.numBeads-1).^(-1.5)))*(1:obj.params.numBeads-1).^(-1.5),...
                  'Color','m',...
                  'LineWidth',3',...
                  'Parent',a,...
                  'DisplayName','expected');
            xlabel(a,'bead distance','FontSize',40);
            ylabel(a,'loop Probability','FontSize',40);
            legend(get(a,'Children'));
        end
    end
end