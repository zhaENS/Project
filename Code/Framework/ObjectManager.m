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
    %TODO: find a uniform interface for the Step of composite and simple objects 
    % TODO: fix minDist for between objects
    
    properties
        numObjects   % number of objects in the class at any time Automatically updated)
        handles         
%         objectList   % keeps the numbers of groups of objects 
%         objectInds   % indices where objects appear in the general list 
        objParams    % parameters for the objects in the simulation 
                
        % Properties af all objects as one 
        particleDist   % pairwise distance between all particles (All)
        connectivity   % connectivity matrix of particles (All)
        fixedParticles % fixed particles (All) 
        curPos         % current particle position (All)
        prevPos        % previous particle position (All)       
        map            % struct('chainNum',cell(1),'inds',cell(1)); % map indices of objects to their chain members
    end
    
    methods 
        function obj = ObjectManager(objectsParams)
            % Class constructor
            obj.numObjects = 0;% initialize with empty objects 
%             obj.objectList = cell(1);% keep the indices of active objects
            obj.objParams  = objectsParams;
            obj.map        = ObjectMapper;
            % register a listener to the countChange event of ObjectMapper
            obj.handles.listener.countChange = addlistener(obj.map,'countChange',@obj.UpdateCount);
        end
        
        function InitializeObjects(obj,domainClass)
            % Initialize the objects 
            % update the general list of properties treating all objects as
            % one 
            cNb = 0;% cumulative number of particles 
            for cIdx = 1:numel(obj.objParams)
                
              % Initizlie Rouse chains 
              obj.handles.chain(cIdx) = Rouse(obj.objParams(cIdx));
              obj.handles.chain(cIdx).SetInitialChainPosition(domainClass);
              inds = (cNb+1):(cNb+obj.objParams(cIdx).numBeads);
              
              % Register the chain as an object in the ObjectMapper class
              obj.map.AddObject(inds)
              
              % Set connectivity 
              obj.connectivity(inds,inds) = obj.handles.chain(cIdx).connectionMap.map;
              
              % set fixed particle num (TODO: consider listing for each
              % object individually)
              obj.fixedParticles          = [obj.fixedParticles, {obj.objParams(cIdx).fixedBeadNum + cNb}];
              
              % update the position list
              obj.prevPos = [obj.prevPos; obj.handles.chain(cIdx).position.prev];
              obj.curPos  = [obj.curPos; obj.handles.chain(cIdx).position.cur];
              
              % Update cumulative object indices 
              cNb = cNb+obj.objParams(cIdx).numBeads;
            end            
              obj.connectivity = logical( obj.connectivity);              
        end
        
        function Merge(obj,objList)
            % Merge all objects in objList into one active object
            % Remove the objects in objectList from the list of active
            % ojects, wrap the objects as one, and add them to the end of
            % the list. 
            obj.map.MergeObjects(objList);
                       
        end
        
        function Break(obj,objNum)
            % A complete splitting of an object to its members 
              obj.map.BreakObject(objNum)
%               obj.numObjects = obj.map.count;
        end
        
        function SplitMember(obj,objNum,memberNum)
            % Split a memeber of the composite object objNum 
            % create a new object for the splitted member
             obj.map.SplitMember(objNum,memberNum);
            
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
                nb = obj.map.GetMemberCount(objList(oIdx),1:obj.map.GetObjectCount(objList(oIdx)));
                c   = curPos(cNb+1:cNb+nb,:);% take only the relevant part for the object
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

                nb = obj.map.GetMemberCount(objList(oIdx),1:obj.map.GetObjectCount(objList(oIdx)));
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
                indList = obj.map.GetObjectMembers(objList(oIdx));% obj.objectList{objList(oIdx)};
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
            inds          = obj.map.GetAllInds(objList);
            connectionMap = obj.connectivity(inds,inds);
            
%             numX = 0;
%             numY = 0;           
%             for o1Idx = 1:numel(objList)
%                  inds1 = obj.map.GetAllInds(objList(o1Idx));%  ObjectMembers(o1Idx);%  obj.objectInds{objList(o1Idx)};
%                 for o2Idx = 1:numel(objList)
%                      inds2 = obj.map.GetAllInds(objList(o2Idx));%obj.objectInds{objList(o2Idx)};
%                       connectionMap((numX+1):(numX+numel(inds1)),(numY+1):(numY+numel(inds2))) = obj.connectivity(inds1,inds2);
%                     numY = numY+numel(inds2);
%                 end
%                  numX = numX+numel(inds1);
%                  numY = 0;
%             end
        end
        
        function fixedParticles = GetFixedParticles(obj,objList)
            % get the list of fixed particles 
            memberList     = obj.map.GetObjectMembers(objList);
            fixedParticles = [obj.fixedParticles{memberList}];
        end
        
        function particleDistance = GetParticleDistance(obj,objList)
            % Get the pairwise distances of particles in objList
            % the objects indices in objList corrospond to objects listed
            % in obj.objectList
%             objList = sort(objList);
            inds = obj.map.GetAllInds(objList);
            particleDistance = obj.particleDist(inds,inds);
            
%             particleDistance = [];
%             numX = 0;
%             numY = 0;
%             for o1Idx = 1:numel(objList)
%                 inds1 = obj.map.GetAllInds(objList(o1Idx));%  obj.objectInds{objList(o1Idx)};
%                 for o2Idx = 1:numel(objList)
%                     % take the indices corresponding to the object
%                     inds2 = obj.map.GetAllInds(objList(o2Idx));%    obj.objectInds{objList(o2Idx)};
%                     particleDistance(numX+1:numX+numel(inds1),numY+1:numY+numel(inds2)) = obj.particleDist(inds1,inds2);
%                     numY = numY+numel(inds2);
%                 end
%                 numX = numX+numel(inds1);
%                 numY = 0;
%             end            
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
                
        function Step(obj,objNum)%TODO: find a uniform interface for composite and simple structures
            % Advance the objects one step in the simulation, apply forces
            % etc. 
            obj.particleDist = ForceManager.GetParticleDistance(obj.curPos);
                         
            for oIdx = 1:numel(objNum)% for each object
                objList = obj.map.GetObjectMembers(oIdx);% members of the objects
               
                if numel(objList)==1
                        % get the indices for the members of the object 
                        beadInds     = obj.map.GetMemberInds(objNum(oIdx),1);
                        beadDistance = obj.particleDist(beadInds,beadInds);
                        % advance each member one step
                        obj.handles.chain(objList(1)).Step(beadDistance)   
                        % update the curPos list 
                        obj.curPos(beadInds,:)  = obj.handles.chain(objList(1)).position.cur;

                else
                   % For composite object made of several sub-objects
                   connectivityMap  = obj.GetConnectivityMapAsOne(objNum(oIdx));
                   particleDistance = obj.GetParticleDistance(objNum(oIdx)); 
                   [~,curMemberPos] = obj.GetMembersPosition(objNum(oIdx));
                   springConst      = obj.GetSpringConstAsOne(objNum(oIdx));
                   minParticleDist  = obj.GetMinParticleDistAsOne(objNum(oIdx));
                   fixedParticleNum = obj.GetFixedParticles(objNum(oIdx));
                   
                   % split parameters to pass to the forceManager
                   par = obj.GetObjectParameters(objNum(oIdx));
                   par = par{1};
                   fp  = [par.forceParams];
                   
                   %TODO: work out fixed bead num for composite objects
                  
                   newPos = ForceManager.ApplyCompositeInternalForces(curMemberPos,particleDistance,connectivityMap,...
                                                         [fp.springForce],[fp.bendingElasticityForce],...
                                                         springConst,[fp.bendingConst],...                                                         
                                                         minParticleDist,fixedParticleNum,0.1);
                                                                                     
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
            objList     = obj.map.GetObjectMembers(objNum);% members of the object             
            cNb         = obj.map.GetMemberCount(objNum,1:obj.map.GetObjectCount(objNum));
            springConst = zeros(cNb);            
            numX        = 0;
            numY        = 0;
            
            for o1Idx = 1:numel(objList)
                 numParticles1 = obj.map.GetMemberCount(objNum,o1Idx);
                for o2Idx = 1:numel(objList)
                    numParticles2 = obj.map.GetMemberCount(objNum,o2Idx);                     
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
        
        function minParticleDist = GetMinParticleDistAsOne(obj,objNum)
            % Get minParticleDist for a composite structure as one big
            % matrix 
            memberList = obj.map.GetObjectMembers(objNum);% members of the object

            cNb             = obj.map.GetMemberCount(objNum,1:obj.map.GetObjectCount(objNum));
            minParticleDist = zeros(cNb);            
            numX = 0;
            numY = 0;
            for o1Idx = 1:numel(memberList)
                 numParticles1 = obj.map.GetMemberCount(objNum,o1Idx);
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
        
        function UpdateCount(obj,sourceObj,varargin)
            % Update number of objects property callback to the listener to
            % the event countChange
            obj.numObjects = sourceObj.count;
        end
        
        function ObjectInteraction(obj)
            % Check for possible interaction between objects and update
            % their data accordingly 
            
            
%             objList= 1:obj.map.count;
% %             %=== test dynamic connectivity ====
%             prob = 0.990;
%             for oIdx = 1:obj.numObjects
%                 r = rand(1);
%                 c = randperm(obj.map.GetMemberCount(oIdx,1:obj.map.GetObjectCount(oIdx)));
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
% %             % ==================================
            
%             % ============ test merging structures ==========
                prob = 0.95;
                r = rand(1);
                o = 1:obj.map.count;% randperm(obj.numObjects);
                obj1 = 1;
               
                if r>0.99
                     if obj.numObjects>1
                          obj2 = 2;% o(end-1);

%                     % the the head of the last object and the tail of the
                      % previous one should be connected 
                      % get indices 
                      
                      indsEnd  = obj.map.GetAllInds(obj1);%obj.objectInds{obj.map.count};
                      indsBend = obj.map.GetAllInds(obj2);%obj.objectInds{obj.map.count-1};
                      obj.connectivity(indsBend(1),indsEnd(end)) = true;
                      obj.connectivity(indsEnd(end),indsBend(1)) = true;
                      
                      obj.Merge([obj1 obj2]);
                      disp('merge')
                     else
               
                     end
                elseif r<(0.01)
                        m = obj.map.GetObjectCount(1:obj.numObjects);
                         % split the biggest one 
                        [numMem, oInd] = max(m);
                        if numMem>1
                        indsEnd   = obj.map.GetMemberInds(oInd,1);%obj.objectInds{obj.map.count};
                        indsBend  = obj.map.GetMemberInds(oInd,2);%obj.objectInds{obj.map.count};
                       obj.connectivity(indsBend(1),indsEnd(end)) = false;
                       obj.connectivity(indsEnd(end),indsBend(1)) = false;
                      
                        obj.SplitMember(oInd,obj.map.GetObjectCount(oInd))
                        
                        end
               disp('split')
%                 elseif r<(1-prob)
%                     % remove connectivity                     
%                      
%                      % disconnect the last two members of the last object 
%                      objCount = obj.map.GetObjectCount(obj2);
%                      if objCount >1
%                         % get the last two members 
%                        indsEnd = obj.map.GetMemberInds(obj2,[(objCount-1), objCount]);
% %                       indsEnd  = obj.map.GetAllInds(obj.map.count);
%                       obj.connectivity(indsEnd(end),indsEnd(1)) = false;
%                       obj.connectivity(indsEnd(1),indsEnd(end)) = false;
%                       % split the last member added
%                       obj.SplitMember(obj2,objCount);
%                        disp('split')
                    
                   
              

            % ============================================
%             encounterDist = 0.1;
            % search for close beads and connect them, update the
            % connectivity accordingly 
%             obj.connectivity = obj.connectivity | obj.particleDist<encounterDist;
          end
        end
    end
end