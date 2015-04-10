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
        tempObjParams = ChainParams;
        particleDist % pairwise distance between all particles
        connectivity % connectivity matrix of all particles
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
%                 [prev{oIdx},cur{oIdx},~,~] = obj.GetObjectsAsOne(oIdx);
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
                c = curPos(cNb+1:cNb+nb,:);% take only the relevant part for the object
                obj.SetCurrentParticlePosition(objList(oIdx),c);
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
                obj.SetPreviousParticlePosition(objList(oIdx),p);
                cNb = cNb+nb;% increase cummulative count 
                                
            end
        end
        
        function SetCurrentParticlePosition(obj,objNum,curPos)
            % objNum is the object number as appears in the objectList. It
            % could be a group of objects. currently only one integer. 
            % Divide the curPos matrix between the members of the
            % object group by their order of appearance in the group 
            
            memberList = obj.objectList{objNum};% members of the object
            cNb        = 0;% cummulative number of object nodes 
            for oIdx = 1:numel(memberList)                                
                oNb = obj.handles.chain(memberList(oIdx)).params.numBeads;
                obj.handles.chain(memberList(oIdx)).position.cur = curPos(cNb+1:cNb+oNb,:);
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
            [~,curPos,obj.connectivity,~] = obj.GetObjectsAsOne(1:obj.numObjects);
            obj.particleDist = ForceManager.GetParticleDistance(curPos);
            for oIdx = 1:numel(objNum)
                objList = obj.objectList{objNum(oIdx)};% members of the objects
                for pIdx = 1:numel(objList)
                    obj.handles.chain(objList(pIdx)).Step% advance each member one step
                end
            end
        end
        
        function [prevPosition, curPosition, connectivityMap,params] = GetObjectsAsOne(obj,objList)% TODO: work out joining chains with different internal forces
            % Arrange the information of the objects into one structure to
            % be passed to Domain and apply global forces on them. will
            % also be used when merging objects 
            cNb = 0; % cummulative number of beads
%             nCb = 0; % number of connnected beads
            for oIdx = 1:numel(objList)
                  indList = obj.objectList{objList(oIdx)};
                  for pIdx = 1:numel(indList)
                      cNb = cNb+ obj.handles.chain(indList(pIdx)).params.numBeads; 
                  end
            end
            
            prevPosition    = zeros(cNb,3);
            curPosition     = zeros(cNb,3);
            connectivityMap = false(cNb);
            cBeads          = []; % cumulative connected beads 
            springConst     = zeros(cNb); % spring constants 
            params          = obj.tempObjParams; 
            fixedBeads      = [];
            cNumBeads       = 0; %cumulative number of beads 
            for oIdx = 1:numel(objList)
                  indList = obj.objectList{objList(oIdx)};
                  for pIdx = 1:numel(indList)
                      numBeads  = obj.handles.chain(indList(pIdx)).params.numBeads;
                      pos       = obj.handles.chain(indList(pIdx)).position;
                      prevPosition(cNumBeads+1:cNumBeads+numBeads,:) = pos.prev;
                      curPosition(cNumBeads+1:cNumBeads+numBeads,:)  = pos.cur;
                      connectivityMap(cNumBeads+1:cNumBeads+numBeads,cNumBeads+1:cNumBeads+numBeads) =...
                          obj.handles.chain(indList(pIdx)).connectionMap.map;
                      connectedBeads = obj.handles.chain(indList(pIdx)).params.connectedBeads;
                      if ~isempty(connectedBeads)
                          cBeads(end+1:end+size(connectedBeads,1),:)         = connectedBeads + cNumBeads;
                      end
                      springConst(cNumBeads+1:cNumBeads+numBeads,cNumBeads+1:cNumBeads+numBeads) = ...
                          obj.handles.chain(indList(pIdx)).params.springConst;
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
            connectivityMap       = logical(connectivityMap);
            
        end
        
        function ObjectInteraction(obj)
            % Check for possible interaction between objects and update
            % their data accordingly 
            encounterDist = 0.1;
            % search for close beads and connect them, update the
            % connectivity accordingly 
            obj.connectivity = obj.connectivity | obj.particleDist<encounterDist;
            
        end
    end
end