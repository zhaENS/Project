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
              obj.handles.chain(cIdx).SetInitialChainPosition(domainClass);
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
        
        function [prev,cur] = GetPosition(obj,objList)
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
        
        function params = GetObjectParameters(obj,objList)
            % get parameters for objects 
            params = cell(1,numel(objList));
            for oIdx = 1:numel(objList)
                % Get all positions for the current group 
                indList = obj.objectList{objList(oIdx)};
                for pIdx = 1:numel(indList)
                  params{oIdx} = [params{oIdx};obj.handles.chain(indList(pIdx)).params];
                end
            end
        end
        
        function DealCurrentPosition(obj,objList,curPos)
            % deal the curPos to the objects and their members 
            cNb  = 0;% cummulative number of beads
            for oIdx = 1:numel(objList)
                % get the number of beads for the object, cut the curPos
                % and send it to SetCurrentPosition with the objNum
                
                nb      = 0;
                indList = obj.objectList{objList(oIdx)};
                for nIdx = 1:numel(indList)
                 nb   = nb+obj.handles.chain(indList).params.numBeads;
                end
                c = curPos(cNb+1:cNb+nb,:);
                obj.SetCurrentParticlePosition(oIdx,c);
                cNb = cNb+nb;% increase cummulative count 
                                
            end
        end
        
        function DealPreviousPosition(obj,objList,prevPos)
            % deal the curPos to the objects and their members 
            cNb  = 0;% cummulative number of beads
            for oIdx = 1:numel(objList)
                % get the number of beads for the object, cut the curPos
                % and send it to SetCurrentPosition with the objNum
                
                nb      = 0;
                indList = obj.objectList{objList(oIdx)};
                for nIdx = 1:numel(indList)
                 nb   = nb+obj.handles.chain(indList).params.numBeads;
                end
                p = prevPos(cNb+1:cNb+nb,:);
                obj.SetPreviousParticlePosition(oIdx,p);
                cNb = cNb+nb;% increase cummulative count 
                                
            end
        end
        
        function SetCurrentParticlePosition(obj,objNum,curPos)
            % objNum is the object number as appears in the objectList. It
            % could be a group of objects. currently an integer. 
            % Divide the curPos matrix between the members of the
            % object group by their order of appearance in the group 
            objList = (obj.objectList{objNum});
            cNb     = 0;% cummulative number of beads 
            for oIdx = 1:numel(objList)                                
                oNb = obj.handles.chain(objList(oIdx)).params.numBeads;
                obj.handles.chain(objList(oIdx)).position.cur = curPos(cNb+1:cNb+oNb,:);
                cNb = cNb+oNb;% increase cummulative count 
            end            
        end
        
        function SetPreviousParticlePosition(obj,objNum,prevPos)
            % objNum is the object number as appears in the objectList. It
            % could be a group of objects. currently an integer. 
            % Divide the curPos matrix between the members of the
            % object group by their order of appearance in the group 
            objList = (obj.objectList{objNum});
            cNb     = 0;% cummulative number of beads 
            for oIdx = 1:numel(objList)                                
                oNb = obj.handles.chain(objList(oIdx)).params.numBeads;
                obj.handles.chain(objList(oIdx)).position.prev = prevPos(cNb+1:cNb+oNb,:);
                cNb = cNb+oNb;% increase cummulative count 
            end
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
        
        function ConnectParticles(obj,objNum,particle1,particle2)
           obj.handles.chain(objNum).ConnectBeads(particle1,particle2);
        end
        
        function DisconnectParticles(obj,objNum,particle1,particle2)
            obj.handles.chain(objNum).DisconnectBeads(particle1,particle2);
        end
                
        function Step(obj,objNum)
            % Advance the objects one step in the simulation , apply forces
            % etc. 
            objList = obj.objectList{objNum};
            for oIdx = 1:numel(objList)
                obj.handles.chain(objList(oIdx)).Step
            end
        end
        
        function [prevPosition, curPosition, connectivityMap,params] = GetObjectsAsOne(obj,objList)% TODO: work out joining chains with different internal forces
            % Arrange the information of the objects into one structure to
            % be passed to Domain and apply global forces on them. will
            % also be used when merging objects 
            prevPosition    = [];
            curPosition     = [];
            connectivityMap = [];
            cBeads          = []; % cumulative connected beads 
            springConst     = []; % spring constants 
            params          = ChainParams;
            fixedBeads      = [];
            cNumBeads       = 0; %cumulative number of beads 
            for oIdx = 1:numel(objList)
                  indList = obj.objectList{objList(oIdx)};
                  for pIdx = 1:numel(indList)
                      numBeads  = obj.handles.chain(indList(pIdx)).params.numBeads;                      
                      prevPosition(end+1:end+numBeads,:) = [obj.handles.chain(indList(pIdx)).position.prev];
                      curPosition(end+1:end+numBeads,:)  = [obj.handles.chain(indList(pIdx)).position.cur];
                      connectivityMap(end+1:end+numBeads,end+1:end+numBeads) = obj.handles.chain(indList(pIdx)).connectionMap.map;
                      connectedBeads = obj.handles.chain(indList(pIdx)).params.connectedBeads;
                      cBeads(end+1:end+size(connectedBeads,1),:)         = connectedBeads + cNumBeads;
                      springConst(end+1:end+numBeads,end+1:end+numBeads) = obj.handles.chain(indList(pIdx)).params.springConst;
                      fBeads    = obj.handles.chain(indList(pIdx)).params.fixedBeadNum;
                      fixedBeads(end+1:end+numel(fBeads)) = fBeads;
                      cNumBeads = cNumBeads+numBeads;
                  end
            end  
            % Adjust parameters
            params.numBeads       = cNumBeads;   
            params.springConst    = springConst;
            params.fixedBeadNum   = fixedBeads;
            params.connectedBeads = cBeads;
            connectivityMap = logical(connectivityMap);
            
        end
        
    end
end