classdef ObjectManager<handle
    % This class manages all objects located inside the domain in the
    % simulation Framework class. 
    % the class is responsible for advancing all objects one step in the
    % simulation, updating their position, keeping track of particle
    % distances, allowing access to information related to a single or
    % multiple objects. and merging several objects into one or breaking
    % objects into several objects 
    properties
        numObjects % number of objects in the class at any time
        handles         
        objectList % keeps the numbers of groups of objects 
        objParams  % parameters for the objects in the simulation 
    end
    
    methods 
        function obj = ObjectManager(objectsParams)
            % Class constructor
            obj.numObjects = 0;% initialize with empty objects 
            obj.objectList = cell(1);% keep the indices of active objects
            obj.objParams  = objectsParams;
        end
        
        function InitializeObjects(obj,domainClass)
            % Initialize the objects 
            for cIdx = 1:numel(obj.objParams)
              % Initizlie Rouse chains 
              obj.handles.chain(cIdx) = Rouse(obj.objParams(cIdx));
              obj.handles.chain(cIdx).SetInitialChainPosition(domainClass)
              obj.numObjects                 = obj.numObjects+1;
              obj.objectList{obj.numObjects} = obj.numObjects;
            end
        end        
        
        function Add(obj,newObject)
            % Add a new initialized object to the list of active objects
            % The newObject should be passed as a class             
            obj.numObjects                     = obj.numObjects+1;
            obj.handles.object(obj.numObjects) = newObject;
            obj.objectList{obj.numObjects}     = obj.numObjects;
        end
        
        function Merge(obj,objList)
            % Merge all objects in objList into one active object
            % Remove the objects in objectList from the list of active
            % ojects, wrap the objects as one, and add them to the end of
            % the list. 
            wrapObjects = obj.objects(objList);
            obj.objects = obj.objects(setdiff(1:obj.numObjects,objList));
            
        end
        
        function Break(obj,obj1,objList)
            % break the obj1 into objects in objList
        end
        
        function [prev,cur]= GetPosition(obj,objList)
            % Get positoins on objects in objList
            % the number of position depends on the grouping given in
            % objectList
            prev = cell(1,numel(objList));
            cur  = cell(1,numel(objList)); 
            for oIdx = 1:numel(objList)
                % Get all positions for the current group 
                indList = obj.objectList{objList(oIdx)};
                for pIdx = 1:numel(indList)
                  prev{oIdx} = [prev{oIdx};obj.handles.chain(indList(pIdx)).position.prev];
                  cur{oIdx} = [cur{oIdx};obj.handles.chain(indList(pIdx)).position.cur];
                end
            end
            
        end
        
        function GetCenterOfMas(obj,objList)
            % Center of mass for objects in objList. 
            % The number of results depends on the groups defined by objectList
            
        end
        
        function [connectionMaps] = GetConnectionMap(obj,objList)
            
            connectionMaps = cell(1,numel(objList));            
            for oIdx = 1:numel(objList)
                % Get all positions for the current group 
                indList = obj.objectList{objList(oIdx)};
                for pIdx = 1:numel(indList)
                  connectionMaps{oIdx} = [connectionMaps{oIdx};obj.handles.chain(indList(pIdx)).connectionMap];
                end
            end 
        end
        
        function Step(obj)
            % Advance the objects one step in the simulation , apply forces
            % etc. 
        end
        
        function UpdatePositon(obj,newPos)
            % update the position of objects by newPos
        end
    end
end