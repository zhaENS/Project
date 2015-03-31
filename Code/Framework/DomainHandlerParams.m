classdef DomainHandlerParams<handle
    % this class holds the parameters of the domain handler
    properties
        dimension@double
        domainShape@char % sphere/cylinder/twoPlates/none
        domainWidth@double
        domainHeight@double
        reflectionType@char % current only option: preserveEnergy
        domainCenter@double % position in space for the center default [0 0 0] 3d
        showDomain@logical  % flag to show or hide the domain 
        parentAxes
    end
    
    methods
        function obj = DomainHandlerParams
            obj.SetDefaultParams
        end
        
        function SetDefaultParams(obj)
            obj.domainShape  = 'sphere';
            obj.domainWidth  = 10;
            obj.domainHeight = 20;
            obj.domainCenter = [0 0 0];
            obj.showDomain   = true;
            obj.reflectionType = 'preserveEnergy';
        end
    end
end
