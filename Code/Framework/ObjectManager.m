classdef ObjectManager<handle
    % This class manages all objects located inside the domain in the
    % simulation Framework class. 
    % the class is responsible for advancing all objects one step in the
    % simulation, updating their position, keeping track of particle
    % distances, allowing access to information related to a single or
    % multiple objects. and merging several objects into one or breaking
    % objects into several objects 
    properties
        numObjects@double % number of objects in the class at any time
        
    end
    
    methods 
        function obj = ObjectManager()
            % Class constructor
        end
        
        function Add(obj,newObject)
            % add the newobject to the list of active objects
        end
        
        function Merge(obj,objList)
            % merge the objects in objList into one active object
        end
        
        function Break(obj,obj1,objList)
            % break the obj1 into objects in objList
        end
        
        function GetPosition(obj,objList)
            % Get positoins on objects in objList
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