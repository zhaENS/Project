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
            params.prob2distMethod = 'fitModel';% options: [fitModel,rouse,composite]
            params.distToAnalyze   = 1;
            params.numDistances    = 1;            
            params.beadsToAnalyze  = 1;
            obj.reconstruction = params;
            
        end
        
        function SetSmootherParams(obj)
            % Signal smoothing parameters
            params.nHoodRad  = 10;         % radius of kernel (full size is [2nHoodRad+1,2nHoodRad+1] in the 2D case)
            params.method     = 'gaussian'; % smoothing method
            params.sigma      = 1;  % std of kernel (or degree)
            params.kernel = []; % smoothing kernel 
            params.kernelAngle = pi/2; % rotation of the kernel 
            
            obj.smoothing      = params;
            
        end
        
        function SetChainParams(obj)
            % Rouse chain parameters
            % Make sure the SimpleRouseParams is in the working path 
            params                = SimpleRouseParams;
            params.dt             = 1e-3;
            params.numSteps       = 1000;
            params.noiseSTD       = 0.005;
            params.b              = sqrt(1.5);
            params.diffusionConst = 1;
            params.numSimulations = 1;
            params.dimension      = 3;
            params.recordPath     = true;% for visualisation 
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