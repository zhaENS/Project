classdef ReconstructionParams<handle
    % Polymer reconstruction parameter class
    properties
        % Class parameters
        data
        reconstruction
        smoothing
        chain
        interpolation
        optimization
    end
    
    methods
        
        function obj = ReconstructionParams()
            
            obj.SetDataParams
            obj.SetReconstructorParams;
            obj.SetSmootherParams;
            obj.SetChainParams;
            obj.SetOptimizationParams;
            obj.SetInterpolationParams;            
        end
        
        function SetDataParams(obj)
            % Data organization parameters
            params.dataFolder     = fullfile(pwd,'ExperimentDataAnalysis');     % default data folder
            params.dataFileName   = 'savedAnalysisTADDAndE';                    % name of dataset
            obj.data = params;
        end
        
        function SetReconstructorParams(obj)
            % Polymer reconstruction params
            params.prob2distMethod = 'rouse';% options: [fitModel,rouse,composite]
            params.distToAnalyze   = 1;
            params.numDistances    = 1;            
            params.beadsToAnalyze  = 1;
            obj.reconstruction = params;
            
        end
        
        function SetSmootherParams(obj)
            % Signal smoothing parameters
            params.nHoodRad  = 20;         % radius of kernel (full size is [2nHoodRad+1,2nHoodRad+1] in the 2D case)
            params.method    = 'Bilateral'; % smoothing method
            params.sigma     = [1.5 1.5];  % std of kernel (or degree) for bilateral sigma=[sigmaf, sigmag]
            smoothKernel     = eye(2*params.nHoodRad+1);
            n                = size(smoothKernel,1);
            for rIdx = 1:n
                smoothKernel = smoothKernel+diag(ones(1,n-rIdx)*1/(rIdx+1),-rIdx)+diag(ones(1,n-rIdx)*1/(rIdx+1),rIdx);
            end
            smoothKernel       = smoothKernel./sum(smoothKernel(:));
            smoothKernel       = fliplr(smoothKernel);
            params.kernel      = smoothKernel; % smoothing kernel 
            params.kernelAngle = pi/2; % rotation of the kernel             
            obj.smoothing      = params;
            
        end
        
        function SetChainParams(obj)
            % Rouse chain parameters these are the parameters for the
            % SimpleRouseFramework class
            % Make sure the SimpleRouseParams is in the working path 
            params                = SimpleRouseParams;
            params.dt             = 1e-3;
            params.numSteps       = 2000;
            params.noiseSTD       = sqrt(2*1*params.dt);
            params.b              = sqrt(1.5);
            params.diffusionConst = 1;
            params.numSimulations = 1;
            params.numRounds      = 30;
            params.analyzeResults = false;
            params.dimension      = 3;
            params.recordPath     = false;% for visualisation 
            params.calculateMSD   = false;
            obj.chain             = params;
        end
        
        function SetOptimizationParams(obj)
            % Data fit and optimization parameters
            
            % Used for the fit optimization of the encounter signal 
            params.model  = fittype('(1/sum(x.^(-beta))).*x.^(-beta)'); 
            params.fitOpt = fitoptions(params.model);
            set(params.fitOpt,'Lower',0,'Upper',1.5,'StartPoint',1,'Robust','off');            
            obj.optimization = params;
        end
        
        function SetInterpolationParams(obj)
            % Data interpolation parameters
            params.zeroInterpolationMethod = 'linear'; % see interp function for options 
            obj.interpolation = params;
        end
        
        function ListParams(obj)% unused
            % list and display parameters by classes            
        end
        
    end
end