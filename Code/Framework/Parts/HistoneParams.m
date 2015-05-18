classdef HistoneParams<handle
    % Histone class parameters 
    properties
        numHistones@double
        forceParams = ForceManagerParams;
    end
    
    methods
        
        function obj = HistoneParams(varargin)
            
            obj.ParseInputParams(varargin)
        end
        
        function ParseInputParams(obj,varargin)
            % parse the input parameters
            p = varargin{:};
            if mod(numel(p),2)~=0
                error('name-value pair must be inserted')
            end            
            for pIdx = 1:numel(p)/2
                obj.(p{2*pIdx-1}) = p{2*pIdx};
            end            
        end
    end
end