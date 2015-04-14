classdef ObjectManager<handle
    % This class manages all objects located inside the domain in the
    % simulation Framework class. 
    % the class is responsible for advancing all objects one step in the
    % simulation, updating their position, keeping track of particle
    % distances, allowing access to information related to a single or
    % multiple objects. and merging several objects into one or breaking
    % objects into several objects 
    
    %TODO: make a mapping from particle number to chain number
    properties
        numObjects % number of objects in the class at any time
        handles         
        objectList % keeps the numbers of groups of objects 
        objParams  % parameters for the objects in the simulation 
        tempObjParams = ChainParams;
        tempObj       = Rouse(ChainParams);
        particleDist % pairwise distance between all particles
        connectivity % connectivity matrix of all particles
        objectInds   % indices where objects appear in the general list 
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
            cNb = 0;% cumulative number of particles 
            for cIdx = 1:numel(obj.objParams)
              % Initizlie Rouse chains 
              obj.handles.chain(cIdx) = Rouse(obj.objParams(cIdx));
              obj.handles.chain(cIdx).SetInitialChainPosition(domainClass);              
              obj.numObjects                 = obj.numObjects+1;
              obj.objectList{obj.numObjects} = obj.numObjects;
              inds = (cNb+1):(cNb+obj.objParams(cIdx).numBeads);
              obj.connectivity(inds,inds) = obj.handles.chain(cIdx).connectionMap.map;
              % update object indices 
              obj.objectInds{cIdx} = inds;
              cNb = cNb+obj.objParams(cIdx).numBeads;
            end
            
            [~,curPos,obj.connectivity,~,~] = obj.GetObjectsAsOne(1:obj.numObjects); 
            obj.particleDist = ForceManager.GetParticleDistance(curPos);% get all distances
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
            
            objGroup       = 1:obj.numObjects;
            objStaying     = setdiff(objGroup,objList);
            toMerge        = (obj.objectList(objList));
            obj.objectList = obj.objectList(objStaying);
            % update indices
            indsToMerge    = cell2mat(obj.objectInds(objList));
            obj.objectInds = obj.objectInds(objStaying);
            
            % group the indices in objList 
            obj.objectList{end+1} = cell2mat(toMerge(:)');
            obj.objectInds{end+1} = indsToMerge;
            obj.numObjects = numel(obj.objectList);

        end
        
        function Split(obj,objNum)
            % A complete splitting of an object to its members 
            objGroup       = 1:obj.numObjects;
            objStaying     = setdiff(objGroup,objNum);
            toSplit        = obj.objectList(objNum);
            indsToSplit    = obj.objectInds(objNum);
            obj.objectInds = obj.objectInds(objStaying);
            obj.objectList = obj.objectList(objStaying);
            
            % group the indices in objList 
            obj.objectList(end+1:end+numel(toSplit)) = toSplit;
            obj.objectInds(end+1:end+numel(toSplit)) = indsToSplit;
            obj.numObjects = numel(obj.objectList);
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
        
        function [prev,cur] = GetMembersPosition(obj,objNum)
            % Get positoins o members of the object objNum

            indList = obj.objectList{objNum};
            prev = cell(1,numel(indList));
            cur  = cell(1,numel(indList)); 
           
                % Get all positions for the current group 
                
                for pIdx = 1:numel(indList)
                  prev{pIdx} = obj.handles.chain(indList(pIdx)).position.prev;
                  cur{pIdx} = obj.handles.chain(indList(pIdx)).position.cur;
                end
           
        end 
        
        function [prevPos,curPos] = GetPositionAsOne(obj,objList)
            % Get current and previous position for all objects in objList
            % as one long list
            prevPos = [];
            curPos  = [];
            for oIdx = 1:numel(objList)
                [p, c]  = obj.GetPosition(objList(oIdx));
                prevPos = [prevPos;p{1}];
                curPos  = [curPos;c{1}];
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
                 nb   = nb+obj.handles.chain(indList(nIdx)).params.numBeads;
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
                 nb   = nb+obj.handles.chain(indList(nIdx)).params.numBeads;
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
            % Get individual connectivity map for each object in objList            
            connectionMaps = cell(1,numel(objList));            
            for oIdx = 1:numel(objList)
                % Get all positions for the current group 
                indList = obj.objectList{objList(oIdx)};
                for pIdx = 1:numel(indList)
                  connectionMaps{oIdx} = [connectionMaps{oIdx};obj.handles.chain(indList(pIdx)).connectionMap];
                end
            end 
        end
        
        function [connectionMap] = GetConnectionMapAsOne(obj,objList)
            % Get the connectivity map of objects in objList in one matrix,
            % displaying the connectivity between these objects.            
            % This function returns the connectivity map only after it was first
            % initialized
            
            numX = 0;
            numY = 0;
            for o1Idx = 1:numel(objList)
                 inds1 = obj.objectInds{objList(o1Idx)};
                for o2Idx = 1:numel(objList)
                     inds2 = obj.objectInds{objList(o2Idx)};
                      connectionMap((numX+1):(numX+numel(inds1)),(numY+1):(numY+numel(inds2))) = obj.connectivity(inds1,inds2);
                    numY = numY+numel(inds2);
                end
                 numX = numX+numel(inds1);
                 numY = 0;
            end
        end
        
        function particleDistance = GetParticleDistance(obj,objList)
            % Get the pairwise distances of particles in objList
            % the objects indices in objList corrospond to objects listed
            % in obj.objectList
            objList = sort(objList);
            particleDistance = [];
            numX = 0;
            numY = 0;
            for o1Idx = 1:numel(objList)
                inds1 = obj.objectInds{objList(o1Idx)};
                for o2Idx = 1:numel(objList)
                    % take the indices corresponding to the object
                    inds2 = obj.objectInds{objList(o2Idx)};
                    particleDistance(numX+1:numX+numel(inds1),numY+1:numY+numel(inds2)) = obj.particleDist(inds1,inds2);
                    numY = numY+numel(inds2);
                end
                numX = numX+numel(inds1);
                numY = 0;
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
            [~,curPos]       = obj.GetPositionAsOne(1:obj.numObjects);
            obj.particleDist = ForceManager.GetParticleDistance(curPos);
            
            cNb = 0; % cummulative number of particles 
            for oIdx = 1:numel(objNum)
                objList = obj.objectList{objNum(oIdx)};% members of the objects
               
                if numel(objList)==1
                    for pIdx = 1:numel(objList)
                        beadInds     = (cNb+1):(cNb+obj.handles.chain(objList(pIdx)).params.numBeads);
                        beadDistance = obj.particleDist(beadInds,beadInds);%TODO: examine what happens if the order of objects changes
                        obj.handles.chain(objList(pIdx)).Step(beadDistance)% advance each member one step
                        cNb  = cNb+obj.handles.chain(objList(pIdx)).params.numBeads;
                    end
                else
                    % for composite object made of several sub-objects
                   connectivityMap  = obj.GetConnectionMapAsOne(oIdx);
                   particleDistance = obj.GetParticleDistance(oIdx); 
                   [~,curPos]       = obj.GetMembersPosition(oIdx);
                   springConst      = obj.GetSpringConstAsOne(oIdx);
                   minParticleDist  = obj.GetMinParticleDistAsOne(oIdx);
                   % split parameters to pass to the forceManager
                   par = obj.GetObjectParameters(oIdx);
                   par = par{1};
                   fp  = [par.forceParams];
                   newPos = ForceManager.ApplyCompositeInternalForces(curPos,particleDistance,connectivityMap,...
                                                         [fp.springForce],[fp.bendingElasticityForce],...
                                                         springConst,[fp.bendingConst],...                                                         
                                                         minParticleDist,[],0.1);%TODO: split parameters
                                                     
                   % Deal the new pos to the object 
                   obj.DealCurrentPosition(oIdx,newPos);
                   obj.DealPreviousPosition(oIdx,newPos);
                end
            end
            % Check for possible interaction between objects
            obj.ObjectInteraction;
        end
        
        function springConst = GetSpringConstAsOne(obj,objNum)% TODO: fix springConst for between objects
            % get the spring constant for a composite structure as one big
            % matrix 
            objList = obj.objectList{objNum};% members of the object
            cNb = 0;
            for oIdx = 1:numel(objList)
                cNb = cNb+obj.handles.chain(objList(oIdx)).params.numBeads;
            end
            springConst = zeros(cNb);
            
            numX = 0;
            numY = 0;
            for o1Idx = 1:numel(objList)
%                  inds1 = obj.objectInds{objList(o1Idx)};
                 numParticles1 = obj.handles.chain(objList(o1Idx)).params.numBeads;
                for o2Idx = 1:numel(objList)
                    numParticles2 = obj.handles.chain(objList(o2Idx)).params.numBeads;
                     
                     if o1Idx==o2Idx
                     springConst((numX+1):(numX+numParticles1),(numY+1):(numY+numParticles2)) = ...
                         obj.handles.chain(objList(o2Idx)).params.springConst;
                     else 
                         % place 1 (Temporary)
                         springConst((numX+1):(numX+numParticles1),(numY+1):(numY+numParticles2)) =1;
                     end
                     numY = numY+numParticles2;
                end
                
                 numX = numX+numParticles1;
                 numY = 0;
            end
            
        end
        
        function minParticleDist = GetMinParticleDistAsOne(obj,objNum)% TODO: fix minDist for between objects
             % get the minParticleDist for a composite structure as one big
            % matrix 
            objList = obj.objectList{objNum};% members of the object
            cNb = 0;
            for oIdx = 1:numel(objList)
                cNb = cNb+obj.handles.chain(objList(oIdx)).params.numBeads;
            end
            minParticleDist = zeros(cNb);
            
            numX = 0;
            numY = 0;
            for o1Idx = 1:numel(objList)
                 numParticles1 = obj.handles.chain(objList(o1Idx)).params.numBeads;
                for o2Idx = 1:numel(objList)
                      numParticles2 = obj.handles.chain(objList(o2Idx)).params.numBeads;
                     if o1Idx==o2Idx
                     minParticleDist((numX+1):(numX+numParticles1),(numY+1):(numY+numParticles2)) = ...
                         obj.handles.chain(objList(o2Idx)).params.minBeadDistance ;
                     else 
                         % place 0 (Temporary)
                         minParticleDist ((numX+1):(numX+numParticles1),(numY+1):(numY+numParticles2)) =0;
                     end
                     numY = numY+numParticles2;
                end
                 numX = numX+numParticles1;
                 numY = 0;
            end 
        end
        
        function [prevPosition, curPosition, connectivityMap,particleDist,params] = GetObjectsAsOne(obj,objList)% TODO: work out joining chains with different internal forces
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
            springForce           = false(numel(objList),1);
            bendingEasticityForce = false(numel(objList),1);
            prevPosition    = zeros(cNb,3);
            curPosition     = zeros(cNb,3);
            connectivityMap = false(cNb);
            cBeads          = []; % cumulative connected beads 
            minBeadDist     = zeros(cNb);
            springConst     = zeros(cNb); % spring constants 
            params          = obj.tempObjParams;            
            fixedBeads      = [];
            cNumBeads       = 0; %cumulative number of beads 
            for oIdx = 1:numel(objList)
                  indList = obj.objectList{objList(oIdx)};
                  for pIdx = 1:numel(indList)
                      numBeads  = obj.handles.chain(indList(pIdx)).params.numBeads;
                      pos       = obj.handles.chain(indList(pIdx)).position;
                      inds      = (cNumBeads+1):(cNumBeads+numBeads);
                      prevPosition(inds,:) = pos.prev;
                      curPosition(inds,:)  = pos.cur;
                      connectivityMap(inds,inds) =...
                          obj.handles.chain(indList(pIdx)).connectionMap.map;
                      connectedBeads = obj.handles.chain(indList(pIdx)).params.connectedBeads;
                      if ~isempty(connectedBeads)
                          cBeads(end+1:end+size(connectedBeads,1),:) = connectedBeads + cNumBeads;
                      end
                      springConst(inds,inds) = ...
                          obj.handles.chain(indList(pIdx)).params.springConst;
                      fBeads    = obj.handles.chain(indList(pIdx)).params.fixedBeadNum;
                      fixedBeads(end+1:end+numel(fBeads)) = fBeads+cNumBeads;
                      minBeadDist(inds,inds) = obj.handles.chain(indList(pIdx)).params.minBeadDistance;
                     
                      cNumBeads = cNumBeads+numBeads;
                  end              
            end  
            % Adjust parameters
            params.numBeads       = cNumBeads;   
            params.springConst    = springConst;
            params.forceParams.springConst = springConst;
            params.fixedBeadNum   = fixedBeads;
            params.forceParams.fixedParticleNum = fixedBeads;
            params.connectedBeads = cBeads;
            params.minBeadDistance= minBeadDist;
            connectivityMap       = logical(connectivityMap);
            particleDist          = [];
            if ~isempty(obj.particleDist)
                particleDist = obj.GetParticleDistance(objList);
            end
        end
        
        function ObjectInteraction(obj)
            % Check for possible interaction between objects and update
            % their data accordingly 
            
            
%             objList= 1:obj.numObjects;
% %             %=== test dynamic connectivity ====
%             prob = 0.990;
%             for oIdx = objList
%                 r = rand(1);
%                 if r>prob
%                     obj.ConnectParticles(oIdx,1,64);
%                 elseif r<(1-prob)
%                     obj.DisconnectParticles(oIdx,1,64)                                              
%                 end
%             end
%             % ==================================
            
            % ============ test merging structures ==========
           prob = 0.99;
                r = rand(1);
                if r>prob
                     if obj.numObjects>1
                      obj.connectivity(1,end) = true;
                      obj.connectivity(end,1) = true;                   
                      obj.Merge([1 obj.numObjects]);
                      disp('merge')
                    end
                elseif r<(1-prob)
%                     obj.connectivity(1,end) = false;
%                     obj.connectivity(end,1) = false;
%                     obj.Split(obj.numObjects);  
%                     disp('split')
                end

            % ============================================
%             encounterDist = 0.1;
            % search for close beads and connect them, update the
            % connectivity accordingly 
%             obj.connectivity = obj.connectivity | obj.particleDist<encounterDist;
            
        end
    end
end