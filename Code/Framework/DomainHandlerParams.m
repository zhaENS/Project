classdef DomainHandlerParams<handle
    % this class holds the parameters of the domain handler
    properties
        dimension@double
        domainShape@char % sphere/cylinder/twoPlates/none
        domainWidth@double
        domainHeight@double
        reflectionType@char % current only option: preserveEnergy
        domainCenter@double % position in space for the center default [0 0 0] 3d
        parentAxes
    end
    
    methods
        function obj = DomainHandlerParams
            obj.SetDefaultParams
        end
        
        function SetDefaultParams(obj)
            obj.domainShape  = 'none';
            obj.domainWidth  = 5;
            obj.domainHeight = 5;
            obj.domainCenter = [0 0 0];
            obj.reflectionType = 'preserveEnergy';
        end
    end
end
