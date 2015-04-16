classdef ObjectManager<handle
    % This class manages all objects located within the domain of the
    % simulation Framework class. 
    % ObjectManager class is responsible for advancing all objects one
    % simulation step,
    % updating objects' positions, keeping track of particle
    % distances, allowing access to information related to a single or
    % multiple objects, merging/splitting objects.
    
    %TODO: make a mapping from particle number to chain number
    %TODO: fix connectivity for composite structures 
    
    properties
        numObjects   % number of objects in the class at any time
        handles         
        objectList   % keeps the numbers of groups of objects 
        objectInds   % indices where objects appear in the general list 
        objParams    % parameters for the objects in the simulation 
        
        tempObjParams = ChainParams;
        tempObj       = Rouse(ChainParams);
        
        % Properties af all objects as one 
        particleDist % pairwise distance between all particles (All)
        connectivity % connectivity matrix of particles (All)
        curPos       % current particle position (All)
        prevPos      % previous particle position (All)       
        map      = ObjectMapper; % struct('chainNum',cell(1),'inds',cell(1)); % map indices of objects to their chain members
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
              inds = (cNb+1):(cNb+obj.objParams(cIdx).numBeads);
              
              % register the chain as an object in the mapper 
              obj.map.AddObject(inds)
              
              % Register the chain as an object 
%               obj.map.chainNum(cIdx) = {cIdx}; % chain numbers 
%               obj.map.inds(cIdx)     = inds;   % chain indices in the general representation
              
%               obj.numObjects                 = obj.numObjects+1;
%               obj.objectList{obj.numObjects} = obj.numObjects;
              obj.connectivity(inds,inds) = obj.handles.chain(cIdx).connectionMap.map;
              
              % set initial connectivity map 
%               obj.connectivity(obj.map.inds{cIdx},obj.map.inds{cIdx}) = obj.handles.chain(cIdx).connectionMap.map;
              obj.prevPos = [obj.prevPos; obj.handles.chain(cIdx).position.prev];
              obj.curPos  = [obj.curPos; obj.handles.chain(cIdx).position.cur];    
              % update object indices 
%               obj.objectInds{cIdx} = obj.map.obj.map.inds{cIdx};% to be removed
              cNb = cNb+obj.objParams(cIdx).numBeads;
            end            
              obj.numObjects = obj.map.count;
%             [obj.prevPos,obj.curPos,obj.connectivity,~,~] = obj.GetObjectsAsOne(1:obj.numObjects); 
        end        
        
        function Add(obj,newObject)% move to objectMapper
            % Add a new initialized object to the list of active objects
            % The newObject should be passed as a class             
            obj.numObjects                     = obj.numObjects+1;
            obj.handles.object(obj.numObjects) = newObject;
            obj.objectList{obj.numObjects}     = obj.numObjects;
            % Modify the map 
          
        end
        
        function Merge(obj,objList)% obsolete- move to objectMapper
            % Merge all objects in objList into one active object
            % Remove the objects in objectList from the list of active
            % ojects, wrap the objects as one, and add them to the end of
            % the list. 
            obj.map.MergeObjects(objList);
            obj.numObjects = obj.map.count;
%             objGroup       = 1:obj.numObjects;
%             objStaying     = setdiff(objGroup,objList);
%             toMerge        = obj.objectList(objList);
%             obj.objectList = obj.objectList(objStaying);
%             % update indices
%             indsToMerge    = cell2mat(obj.objectInds(objList));
%             
%             % remove the merged obj from the list 
%             obj.objectInds = obj.objectInds(objStaying);
%             
%             % group the indices in objList 
%             obj.objectList{end+1} = cell2mat(toMerge(:)');
%             obj.objectInds{end+1} = indsToMerge;
%             obj.numObjects = numel(obj.objectList);
                       
        end
        
        function Split(obj,objNum)% obsolete move to objectMapper
            % A complete splitting of an object to its members 
%             objGroup       = 1:obj.numObjects;
%             objStaying     = setdiff(objGroup,objNum);
%             toSplit        = obj.objectList(objNum);
%             indsToSplit    = obj.objectInds(objNum);
%             obj.objectInds = obj.objectInds(objStaying);
%             obj.objectList = obj.objectList(objStaying);
%             
%             % group the indices in objList 
%             obj.objectList(end+1:end+numel(toSplit)) = toSplit;
%             obj.objectInds(end+1:end+numel(toSplit)) = indsToSplit;
%             obj.numObjects = numel(obj.objectList);
              obj.map.SplitObject(objNum)
              obj.numObjects = obj.map.count;
        end
        
        function SplitMember(obj,objNum,memberNum)% obsolete move to objectMapper
            % Split a memeber of the composite object objNum 
            % create a new object for the splitted member
            objGroup = 1:obj.numObjects;
            
        end
        
        function [prev,cur] = GetPosition(obj,objList)
            % Get position of objects in objList
            % the number of position depends on the grouping given in
            % objectList
            prev = cell(1,numel(objList));
            cur  = cell(1,numel(objList)); 
            for oIdx = 1:numel(objList)
                % Get all positions for the current group 
                inds       = obj.map.GetAllInds(objList(oIdx));% obj.map.object(objList(oIdx)).inds.allInds; %objectInds{oIdx};
                prev{oIdx} = obj.prevPos(inds,:);
                cur{oIdx}  = obj.curPos(inds,:);
            end            
        end        
        
        function [prev,cur] = GetMembersPosition(obj,objNum)
            % Get positoins of members of the object objNum            
            memberList = obj.map.GetObjectMembers(objNum);%[obj.map.object(objNum).members];% obj.objectList{objNum};% members indices
            prev    = cell(1,numel(memberList));
            cur     = cell(1,numel(memberList));
            
            % Get all positions for the current group            
            for pIdx = 1:numel(memberList)
                prev{pIdx} = obj.handles.chain(memberList(pIdx)).position.prev;
                cur{pIdx}  = obj.handles.chain(memberList(pIdx)).position.cur;
            end
           
        end 
        
        function [prevPos,curPos] = GetPositionAsOne(obj,objList)
            % Get current and previous position for all objects in objList
            % as one long list
            inds    = obj.map.GetAllInds(objList);% [obj.objectInds{objList}];
            prevPos = obj.prevPos(inds,:);
            curPos  = obj.curPos(inds,:);
            
        end
        
        function params = GetObjectParameters(obj,objList)
            % get parameters for objects 
            params = cell(1,numel(objList));
            for oIdx = 1:numel(objList)
                % Get all positions for the current group 
                memberList = obj.map.GetObjectMembers(objList(oIdx)); %object(objList(oIdx)).members; % obj.objectList{objList(oIdx)};% members
                for pIdx = 1:numel(memberList)
                  params{oIdx} = [params{oIdx};obj.handles.chain(memberList(pIdx)).params];
                end
            end
        end
        
        function DealCurrentPosition(obj,objList,curPos)
            % deal the curPos to the objects and their members 
            % curPos is an Nxdim  matrix of positions to be dealt to the
            % members of the objects in objList 
            % the function deals the position according to the order of
            % appearance in the curPos and the object list 
            
            cNb  = 0;% cummulative number of beads
            for oIdx = 1:numel(objList)
                % get the number of beads for the object, cut the curPos
                % and send it to SetCurrentPosition with the objNum
%                 nb         = 0;
%                 memberList = obj.map.GetMembers(objList(oIdx));% obj.objectList{objList(oIdx)};% members 
%                 
%                 for nIdx = 1:numel(memberList)
%                  nb   = nb+obj.handles.chain(memberList(nIdx)).params.numBeads;
%                 end
                
                nb = obj.map.GetObjectCount(objList(oIdx));% ,numel(obj.map.object(oIdx).inds.allInds);
                c   = curPos(cNb+1:cNb+nb,:);% take only the relevant part for the object
                obj.SetCurrentParticlePosition(objList(oIdx),c);
                
%                 % update the general list 
%                 inds = obj.map.GetAllInds(objList(oIdx));% obj.objectInds{objList(oIdx)};
%                 obj.curPos(inds,:) = c;
                cNb = cNb+nb;% increase cummulative count                                 
            end
        end
        
        function DealPreviousPosition(obj,objList,prevPos)
            % deal the curPos to the objects and their members 
            cNb  = 0;% cummulative number of beads
            for oIdx = 1:numel(objList)
                % get the number of beads for the object, cut the curPos
                % and send it to SetCurrentPosition with the objNum                
%                 nb      = 0;
%                 memberList = obj.map.GetObjectMembers(objList(oIdx));% objectList{objList(oIdx)};
%                 for nIdx = 1:numel(indList)
%                  nb   = nb+obj.handles.chain(indList(nIdx)).params.numBeads;
%                 end
                nb = obj.map.GetObjectCount(objList(oIdx));
                p = prevPos(cNb+1:cNb+nb,:);
                obj.SetPreviousParticlePosition(objList(oIdx),p);
                % update the general list 
%                 inds = obj.map.GetAllInds(objList(oIdx));
% %                 inds = obj.objectInds{objList(oIdx)};
%                 obj.prevPos(inds,:) = p;
                cNb = cNb+nb;% increase cummulative count                                 
            end
        end
        
        function SetCurrentParticlePosition(obj,objNum,curPos)
            % objNum is the object number as appears in the objectList. It
            % could be a group of objects. currently only one integer. 
            % Divide the curPos matrix between the members of the
            % object group by their order of appearance in the group 
            
            memberList = obj.map.GetObjectMembers(objNum);%    obj.objectList{objNum};% members of the object
            cNb        = 0;% cummulative number of object nodes 
            for oIdx = 1:numel(memberList)                                
%                 oNb = obj.handles.chain(memberList(oIdx)).params.numBeads;
                oNb = obj.map.GetMemberCount(objNum,oIdx);
                obj.handles.chain(memberList(oIdx)).position.cur = curPos(cNb+1:cNb+oNb,:);
                cNb = cNb+oNb;% increase cummulative count 
            end            
            
            % update the general list 
            inds = obj.map.GetAllInds(objNum);
            obj.curPos(inds,:) = curPos;
        end
        
        function SetPreviousParticlePosition(obj,objNum,prevPos)
            % objNum is the object number as appears in the objectList. It
            % could be a group of objects. currently an integer. 
            % Divide the curPos matrix between the members of the
            % object group by their order of appearance in the group 
            memberList = obj.map.GetObjectMembers(objNum);% (obj.objectList{objNum});
            cNb        = 0;% cummulative number of beads 
            for oIdx = 1:numel(memberList)                                
%                 oNb = obj.handles.chain(objList(oIdx)).params.numBeads;
                oNb = obj.map.GetMemberCount(objNum,oIdx);
                obj.handles.chain(memberList(oIdx)).position.prev = prevPos(cNb+1:cNb+oNb,:);
                cNb = cNb+oNb;% increase cummulative count 
            end
            
            % Update the general list 
            inds = obj.map.GetAllInds(objNum);
            obj.prevPos(inds,:) = prevPos;
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
        
        function [connectionMap] = GetConnectivityMapAsOne(obj,objList)
            % Get the connectivity map of objects in objList in one matrix,
            % displaying the connectivity between these objects.            
            % This function returns the connectivity map only after it was first
            % initialized
            
            numX = 0;
            numY = 0;           
            for o1Idx = 1:numel(objList)
                 inds1 = obj.map.GetAllInds(objList(o1Idx));%  ObjectMembers(o1Idx);%  obj.objectInds{objList(o1Idx)};
                for o2Idx = 1:numel(objList)
                     inds2 = obj.map.GetAllInds(objList(o2Idx));%obj.objectInds{objList(o2Idx)};
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
%             objList = sort(objList);
            particleDistance = [];
            numX = 0;
            numY = 0;
            for o1Idx = 1:numel(objList)
                inds1 = obj.map.GetAllInds(objList(o1Idx));%  obj.objectInds{objList(o1Idx)};
                for o2Idx = 1:numel(objList)
                    % take the indices corresponding to the object
                    inds2 = obj.map.GetAllInds(objList(o2Idx));%    obj.objectInds{objList(o2Idx)};
                    particleDistance(numX+1:numX+numel(inds1),numY+1:numY+numel(inds2)) = obj.particleDist(inds1,inds2);
                    numY = numY+numel(inds2);
                end
                numX = numX+numel(inds1);
                numY = 0;
            end            
        end
        
        function ConnectParticles(obj,objNum,particle1,particle2)
            % Connect two particles
            % particle1 and particle2 are number of particles in the object            
            objList = obj.map.GetObjectMembers(objNum); %obj.objectList{objNum};% chain members of the object
            
%             % get the chain 
%             map.particleChain(particle1)
%             map.particleChain(particle2)
            
            for pIdx = 1:numel(objList)
               % if the particles belong to the same chain
              obj.handles.chain(objList(pIdx)).ConnectBeads(particle1,particle2);             
                % for composite object change the connectivity             
            end
            
        end
        
        function DisconnectParticles(obj,objNum,particle1,particle2)%TODO: apply to composite objects as well
            % Disconnect pairs of particles in an object objNum
            % the particle indices should be relative to the ones in objNum
            for pIdx = 1:size(particle1,1)
              obj.handles.chain(objNum).DisconnectBeads(particle1(pIdx),particle2(pIdx));
              % Update the general connectivity map 
              
            end
        end
                
        function Step(obj,objNum)
            % Advance the objects one step in the simulation, apply forces
            % etc. 
            obj.particleDist = ForceManager.GetParticleDistance(obj.curPos);
            
            for oIdx = 1:numel(objNum)% for each object
                objList = obj.map.GetObjectMembers(oIdx);% obj.objectList{objNum(oIdx)};% members of the objects
               
                if numel(objList)==1
                        beadInds     = obj.map.GetMemberInds(objNum(oIdx),1);% obj.objectInds{objNum(oIdx)};% get the indices of the object 
                        beadDistance = obj.particleDist(beadInds,beadInds);%TODO: examine what happens if the order of objects changes
                        obj.handles.chain(objList(1)).Step(beadDistance)% advance each member one step
                        % update the curPos list 
                        obj.curPos(beadInds,:)  = obj.handles.chain(objList(1)).position.cur;

                else
                    % For composite object made of several sub-objects
                   connectivityMap  = obj.GetConnectivityMapAsOne(objNum(oIdx));
                   particleDistance = obj.GetParticleDistance(objNum(oIdx)); 
                   [~,curMemberPos] = obj.GetMembersPosition(objNum(oIdx));
                   springConst      = obj.GetSpringConstAsOne(objNum(oIdx));
                   minParticleDist  = obj.GetMinParticleDistAsOne(objNum(oIdx));
                   
                   % split parameters to pass to the forceManager
                   par = obj.GetObjectParameters(objNum(oIdx));
                   par = par{1};
                   fp  = [par.forceParams];
                   
                   %TODO: work out fixed bead num for composite objects
                   fixedBeadNum = [];
                   newPos = ForceManager.ApplyCompositeInternalForces(curMemberPos,particleDistance,connectivityMap,...
                                                         [fp.springForce],[fp.bendingElasticityForce],...
                                                         springConst,[fp.bendingConst],...                                                         
                                                         minParticleDist,fixedBeadNum,0.1);
                                                                                     
                % Deal the new pos to the object 
                   obj.DealCurrentPosition(oIdx,newPos);
                % obj.DealPreviousPosition(oIdx,newPos);
                 end
            end
            
            % Check for possible interaction between objects
            obj.ObjectInteraction;
        end
        
        function springConst = GetSpringConstAsOne(obj,objNum)% TODO: fix springConst for between objects
            % get the spring constant for a composite structure as one big
            % matrix 
            objList = obj.map.GetObjectMembers(objNum);% obj.objectList{objNum};% members of the object
%             cNb = 0;
%             for oIdx = 1:numel(objList)
%                 cNb = cNb+obj.handles.chain(objList(oIdx)).params.numBeads;
%             end
%             
            cNb = numel(obj.map.object(objNum).inds.allInds);
            springConst = zeros(cNb);
            
            numX = 0;
            numY = 0;
            for o1Idx = 1:numel(objList)
%                  inds1 = obj.objectInds{objList(o1Idx)};
                 numParticles1 = obj.map.GetMemberCount(objNum,o1Idx);%numel(obj.map.object(objNum).inds.memberInds(o1Idx)); %obj.handles.chain(objList(o1Idx)).params.numBeads;
                for o2Idx = 1:numel(objList)
                    numParticles2 = obj.map.GetMemberCount(objNum,o2Idx); % numel(obj.map.object(objNum).inds.memberInds(o2Idx));% obj.handles.chain(objList(o2Idx)).params.numBeads;
                     
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
            memberList = obj.map.GetObjectMembers(objNum);% obj.objectList{objNum};% members of the object
%             cNb = 0;
%             for oIdx = 1:numel(objList)
%                 cNb = cNb+obj.handles.chain(objList(oIdx)).params.numBeads;
%             end
            cNb             = obj.map.GetObjectCount(objNum);%numel(obj.map.GetAllInds(objNum));
            minParticleDist = zeros(cNb);            
            numX = 0;
            numY = 0;
            for o1Idx = 1:numel(memberList)
                 numParticles1 = obj.map.GetMemberCount(objNum,o1Idx);% obj.handles.chain(objList(o1Idx)).params.numBeads;
                for o2Idx = 1:numel(memberList)
                      numParticles2 = obj.map.GetMemberCount(objNum,o2Idx);%obj.handles.chain(objList(o2Idx)).params.numBeads;
                     if o1Idx==o2Idx
                     minParticleDist((numX+1):(numX+numParticles1),(numY+1):(numY+numParticles2)) = ...
                         obj.handles.chain(memberList(o2Idx)).params.minBeadDistance ;
                     else 
                         % place 0 (Temporary)
                         minParticleDist ((numX+1):(numX+numParticles1),(numY+1):(numY+numParticles2)) = 0;
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
%                 c = randperm(32);
%                 if r>prob
% %                     cBeads = obj.handles.chain(oIdx).params.connectedBeads;
% %                     if ~isempty(cBeads)
% %                     obj.DisconnectParticles(oIdx,cBeads(:,1),cBeads(:,2))                     
% %                     end
%                     obj.ConnectParticles(oIdx,c(1),c(2));
%                 elseif r<(0.01)
%                     cBeads = obj.handles.chain(oIdx).params.connectedBeads;
%                     if ~isempty(cBeads)
%                     obj.DisconnectParticles(oIdx,cBeads(:,1),cBeads(:,2));
%                     end
%                 end
%             end
%             % ==================================
            
%             % ============ test merging structures ==========
           prob = 0.95;
                r = rand(1);
                if r>prob
                     if obj.map.count>1
                      % choose 2 random particles 
%                       c = randperm(obj.numObjects);
%                       % find the object indices matching the connected
%                       % beads 
%                       obj1 = ceil(c(1)/32);
%                       obj2 = ceil(c(2)/32);
% %                       if c(1)~=c(2)
%                           % choose 2 random beads 
%                           b = randperm(32);
%                           bead1 = (c(1)-1)*32+b(1);
%                           bead2 = (c(2)-1)*32+b(2);
% %                           if obj1==obj2
% %                               obj.handles.chain(obj1).ConnectBeads(
% %                           end
%                     % the the head of the last object and the tail of the
                      % previous one should be connected 
                      % get indices 
                      if obj.map.count~=1
                      indsEnd  = obj.map.GetAllInds(obj.map.count);%obj.objectInds{obj.numObjects};
                      indsBend = obj.map.GetAllInds(obj.map.count-1);%obj.objectInds{obj.numObjects-1};
                      obj.connectivity(indsBend(end),indsEnd(1)) = true;
                      obj.connectivity(indsEnd(1),indsBend(end)) = true;
                      
                      obj.Merge([obj.numObjects-1 obj.numObjects]);
                      disp('merge')
                      end
                    end
                elseif r<(1-prob)
                    
%                     obj.connectivity(1,3*32) = false;
%                     obj.connectivity(96,1) = false;
                    
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