classdef ObjectManager<handle
    % This class manages all objects located within the domain of the
    % simulation Framework class. 
    % ObjectManager class is responsible for advancing all objects one
    % simulation step,
    % updating objects' positions, keeping track of particle
    % distances, allowing access to information related to a single or
    % multiple objects, merging/splitting objects.
    
    %TODO: find a uniform interface for the Step of composite and simple objects 
    
    properties (SetObservable)
        numObjects   % number of objects in the class at any time Automatically updated)
        handles         
%         objectList   % keeps the numbers of groups of objects 
%         objectInds   % indices where objects appear in the general list 
        objParams    % parameters for the objects in the simulation 
                
        % Properties af all objects as one 
        particleDist   % pairwise distance between all particles (All)
        connectivity   % connectivity matrix of particles (All)
        fixedParticles % fixed particles (All) 
        stickyParticles     %sticky particles (All)
        particlesOnBoundary %particles on boundary (All)
        connectedStickyBeads %the num of sticky beads that have already connected
        domainInds     % the indices of domain the particle belong to 
        curPos         % current particle position (All)
        prevPos        % previous particle position (All)       
        map            % struct('chainNum',cell(1),'inds',cell(1)); % map indices of objects to their chain members
    end
    
    events
        curPosChange       % event for changing the current position
        prevPosChange      % event for changing the previous position 
        connectivityChange % event for changing connectivity 
    end
    
    methods 
        function obj = ObjectManager(objectsParams)
            % Class constructor
            obj.numObjects = 0;% initialize with empty objects 
            obj.objParams  = objectsParams;
            obj.map        = ObjectMapper;
            
            % Register a listener to the countChange event of ObjectMapper
            obj.handles.listener.countChange        = addlistener(obj.map,'countChange',@obj.UpdateCount);
%             obj.handles.listener.connectivityChange = addlistener(obj,'connectivity','PostSet',@obj.ConnectivityChangeListenerCallback);
            
        end
        
        function InitializeObjects(obj,domainClass)
            % Initialize the objects 
            % update the general list of properties treating all objects as
            % one 
            cNb = 0;% cumulative number of particles 
            for cIdx = 1:numel(obj.objParams)
                
              % Initizlie Rouse chains 
              inds = (cNb+1):(cNb+obj.objParams(cIdx).numBeads);
              obj.handles.chain(cIdx) = Rouse(obj.objParams(cIdx),inds,obj);
              obj.handles.chain(cIdx).SetInitialChainPosition(domainClass);
                          
              % Register the chain as an object in the ObjectMapper class
              obj.map.AddObject(inds)
              
              % Set connectivity 
              obj.connectivity(inds,inds) = obj.handles.chain(cIdx).connectionMap.map;
              
              % domain indices 
              obj.domainInds(inds,1) = obj.objParams(cIdx).initializeInDomain;
              
              % set fixed particle num (TODO: consider listing for each
              % object individually)
              obj.fixedParticles          = [obj.fixedParticles, {obj.objParams(cIdx).fixedBeadNum + cNb}];
              
              %set sticky particle num
              obj.stickyParticles         = [obj.stickyParticles, {obj.objParams(cIdx).stickyBeads + cNb}];
              
              % update the position list
              obj.prevPos = [obj.prevPos; obj.handles.chain(cIdx).position.prev];
              obj.curPos  = [obj.curPos; obj.handles.chain(cIdx).position.cur];
         
              % Update cumulative object indices 
              cNb = cNb+obj.objParams(cIdx).numBeads;
            end 
              obj.connectivity = (logical(obj.connectivity));
          
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
        end
        
        function SplitMember(obj,objNum,memberNum)
            % Split a memeber of the composite object objNum 
            % create a new object for the splitted member
             obj.map.SplitMember(objNum,memberNum);
            
        end
        
        function [prev,cur] = GetPosition(obj,objList)%obsolete
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
            memberList = (obj.map.GetObjectMembers(objNum));% members indices
            prev    = cell(1,numel(memberList));
            cur     = cell(1,numel(memberList));            
            
            % Get all positions for the current group            
            for pIdx = 1:numel(memberList)
                % get member indices in the general map
                inds       = obj.map.GetMemberInds(objNum,pIdx);
                prev{pIdx} = obj.prevPos(inds,:);
                cur{pIdx}  = obj.curPos(inds,:);

            end
           
        end 
        
        function [prevPos,curPos] = GetPositionAsOne(obj,objList)
            % Get current and previous position for all objects in objList
            % as one long list
            objList = (objList);
            inds    = obj.map.GetAllInds(objList);% [obj.objectInds{objList}];
            prevPos = obj.prevPos(inds,:);
            curPos  = obj.curPos(inds,:);            
        end
        
        function params = GetObjectParameters(obj,objList)
            % get parameters for objects 
            params = cell(1,numel(objList));
            for oIdx = 1:numel(objList)
                % Get all positions for the current group 
                memberList = obj.map.GetObjectMembers(objList(oIdx)); % members
                for pIdx = 1:numel(memberList)
                  params{oIdx} = [params{oIdx};obj.handles.chain(memberList(pIdx)).params];
                end
            end
        end
        
        function DealCurrentPosition(obj,objList,curPos)
            % Deal the curPos to the objects and their members 
            % curPos is an Nxdim  matrix of positions to be dealt to the
            % members of the objects in objList 
            % the function deals the position according to the order of
            % appearance in the curPos and the object list             
            
            inds = (obj.map.GetAllInds(objList));
            obj.curPos(inds,:) = curPos;% update the position
            notify(obj,'curPosChange'); % notify all registered listeners

        end
        
        function DealPreviousPosition(obj,objList,prevPos)
          %   Update the general list 
          %   deal the prevPos to all the object members

             inds = (obj.map.GetAllInds(objList));
             obj.prevPos(inds,:)= prevPos;
             notify(obj,'prevPosChange'); % notify all registered listeners
            
        end
        
        function SetCurrentParticlePosition(obj,objNum,curPos)%unused
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
        
        function SetPreviousParticlePosition(obj,objNum,prevPos)%unused
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
            objList       = (objList);
            inds          = obj.map.GetAllInds(objList);
            connectionMap = obj.connectivity(inds,inds);

        end
        
       function fixedParticles  = GetFixedParticles(obj,objList)
            % get the list of fixed particles 
            memberList     = obj.map.GetObjectMembers(objList);
            fixedParticles = [obj.objParams(memberList).fixedBeadNum];
            if ~isempty(fixedParticles)
                for mIdx = 2:numel(memberList)
                    fixedParticles(mIdx) = fixedParticles(mIdx)+obj.objParams(memberList(mIdx-1)).numBeads;
                end
            end
            % fixedParticles = [obj.fixedParticles{memberList}];
        end
        
        function stickyParticles = GetStickyParticles(obj,objList)
            % get the list of sticky particles 
            memberList      = obj.map.GetObjectMembers(objList);
            stickyParticles = [obj.objParams(memberList).stickyBeads];
            if ~isempty(stickyParticles)
                for mIdx = 2:numel(memberList)
                    stickyParticles(mIdx) = stickyParticles(mIdx)+obj.objParams(memberList(mIdx-1)).numBeads;
                end
            end
        end
        
        function particleOnBoundary = GetParticlesOnBoundary(obj,objList)
            % get the list of particles on the boundary of (any) domain
            memberList      = obj.map.GetObjectMembers(objList);
            particleOnBoundary = [obj.objParams(memberList).beadsOnBoundary];
            if ~isempty(particleOnBoundary) && (numel(memberList)>1)
                for lIdx = 1:numel(memberList)-1
                    pSize1(lIdx)  = numel(obj.objParams(memberList(lIdx)).beadsOnBoundary);
                    pSize2(lIdx) = numel(obj.objParams(memberList(lIdx+1)).beadsOnBoundary);
                    mIdx         = sum(pSize1(1:lIdx))+1 : sum(pSize2(lIdx))+sum(pSize1(1:lIdx));
                    particleOnBoundary(mIdx) = particleOnBoundary(mIdx)+sum([obj.objParams(memberList(1:lIdx)).numBeads]);
                    
                end
            end
        
        end
            
        function particleDistance = GetParticleDistance(obj,objList)
            % Get the pairwise distances of particles in objList
            % the objects indices in objList corrospond to objects listed
            % in obj.objectList
%             objList = sort(objList);
            inds = obj.map.GetAllInds(objList);
            particleDistance = obj.particleDist(inds,inds);
         
        end
        
        function connectedParticles = GetObjectConnectedParticles(obj,objNum,connectivityState)
            % Get indices of connected particles (locally) for objNum
            % connectivityState is 'full' or 'offDiagonals'. full gives the
            % full map including nearesr neighbor connectivity, whereas
            % 'offDiagonals' gives only off super and sub diagonal
            % connectivity. For a composite structure in which the ends con
            % also be connected the function 
            
            if nargin==2
                connectivityState = 'full';
            end
            
            inds = obj.map.GetAllInds(objNum);
            if strcmpi(connectivityState,'full')
                   cMat = obj.connectivity(inds,inds);
                   [connectedParticles(:,1),connectedParticles(:,2)]= obj.ListConnectivity(cMat);
            
            elseif strcmpi(connectivityState,'offDiagonals')
                % first find off diagonal 
                nnConnectivity = ~diag(true(1,numel(inds)-1),1);
%                 % remove the boundaries between objects from the diagonal 
%                 c = obj.map.object(objNum).inds.memberCount;% GetObjectCount(objNum);
%                 for cIdx = 1:numel(c)
%                     % Remove the last line of the nnConnectivity (no
%                     % neighor with higher ind is connected to the 'right' of
%                     % the last ind)
%                     nnConnectivity(cumsum(c(1:cIdx)),:)= false;
%                 end                
                cMat = obj.connectivity(inds,inds)& nnConnectivity;
                [connectedParticles(:,1),connectedParticles(:,2)] = obj.ListConnectivity(cMat);                
            end
            
            if ~isempty(connectedParticles)
                connectedParticles(:,1) = obj.map.object(objNum).inds.allInds(connectedParticles(:,1));
                connectedParticles(:,2) = obj.map.object(objNum).inds.allInds(connectedParticles(:,2));
            end
        end
        
        function connectionMap = GetMembersConnectivityMap(obj,memberList)
            % Get the connectivity map between two members. the member
            % numbers are the ones appearing in the global member list 
           
            % Find the objects the members belong to
            memberList = (memberList);
            objNum = obj.map.GetObjectFromMember(memberList);
            
            % Find the indices of the members in these objects 
            inds = [];
            for oIdx = 1:numel(objNum)
                membNum = obj.map.GetMemberNumInObject(objNum(oIdx),memberList(oIdx));
                inds    = [inds, obj.map.GetMemberInds(objNum(oIdx),membNum)];
            end
             connectionMap = obj.connectivity(inds,inds);
                        
        end
        
        function connectedParticles = GetMembersConnectedParticles(obj,memberList,connectivityState)
                % List the connected particles 
                  cMat = obj.GetMembersConnectivityMap(memberList);
            if nargin==2
              connectivityState = 'full'; % include nearest neighbors 
            end

            if strcmpi(connectivityState,'offDiagonals')
                  cMat = cMat & ~diag(true(1,size(cMat,1)-1),1);
            end
                  [connectedParticles(:,1), connectedParticles(:,2)] = obj.ListConnectivity(cMat);
        end
        
        function objDistance = GetObjectDistance(obj,objList)
            % Get the pairwise distance between all members of the
            % specified objects in objlist 9 by the order of their
            % appearance
            inds = obj.map.GetAllInds(objList);
            objDistance = obj.particleDist(inds,inds);
            
        end
        
        function memberDist = GetMemberDistance(obj,memberList)%unfinished
            %get the pairwise distances between members in memberList
            memberDist = [];
            objs = obj.map.GetObjectFromMember(memberList);
            inds = obj.map.GetAllInds();
        end        
        
        
         function AttachToBoundary(obj,encounterDistance)
            %allow beads in domain to attach to boundary if distance is
            %below encounter distance;
            AttachProbability = obj.objParams.probAttachToBoundary;
            for oIdx=1:obj.numObjects
                beadsOnBoundary    = obj.GetParticlesOnBoundary(oIdx);
                [prevPosition,curPosition]    = obj.GetPosition(oIdx);
                prevPosition           = prevPosition{1};
                curPosition        = curPosition{1};
                D                  = pdist2mex(curPosition',curPosition','euc',[],[],[],[]);
                [row,col]          = find(D(beadsOnBoundary,:)<encounterDistance & D(beadsOnBoundary,:)>0);
                if ~isempty(col)
                    for cIdx = 1:numel(beadsOnBoundary)
                        col(find(beadsOnBoundary(cIdx)==col))=[];
                    end
                    r   = rand(numel(col),1);
                    attachIndx = r<AttachProbability;
                    for rIdx = 1:numel(col)
                        if attachIndx(rIdx)
                            sprintf('Attachment to the boundary! beads %d of chain %d',col(rIdx),oIdx)
                            %update the position such that they are on the
                            %boundary;
                            [obj.handles.chain(oIdx).params.beadsOnBoundary] = ...
                            [obj.handles.chain(oIdx).params.beadsOnBoundary col(rIdx)];
                            
                            obj.handles.chain(oIdx).params.beadsOnBoundary = ...
                            sort(obj.handles.chain(oIdx).params.beadsOnBoundary);
                            curPosition(col(rIdx),:)  = curPosition(beadsOnBoundary(row(rIdx)),:);
                            prevPosition(col(rIdx),:) = prevPosition(beadsOnBoundary(row(rIdx)),:);
                            obj.DealCurrentPosition(oIdx,curPosition);                              
                            obj.DealPreviousPosition(oIdx,prevPosition);
                        end
                    end
                end
            end         
        end
        
        
        function ConnectStickyParticles(obj,stickyDistance)
                  %calculate each pair of beads in the list stickyBeads to see
                  %if their distances are smaller than encounterDist;
                  %return objNum to see if they are in the same object or
                  %not;
                  stickyBeads = [obj.stickyParticles{:}];
                  A          = obj.curPos(stickyBeads,:);
                  D          = pdist2mex(A',A','euc',[],[],[],[]);
                  [row ,col] = find(tril(D)< stickyDistance & tril(D)>0);
                  % get the probability for two sticky beads to attach
                  stickyBeadConnectProb    =  obj.objParams.probAttachToStickyBeads;% the attachment prob.
                  stickyBeadDisconnectProb = 1-stickyBeadConnectProb;
                  r = rand(numel(row),1);%random values to test if the beads should be connected or disconnected;
                  connectInd = r<stickyBeadConnectProb;
                  if ~isempty(row)
                       for rIdx = 1:numel(row)
                          % Test if two sticky beads should be connected
                          if connectInd(rIdx)
                                sprintf('sticky')
                                sprintf('%d and %d',stickyBeads(row(rIdx)),stickyBeads(row(rIdx)))          
                              obj.ConnectParticles(stickyBeads(row(rIdx)),stickyBeads(col(rIdx)));
                           end
                      end
                  end
                  
                  for oIdx=1:obj.numObjects
                  %Get the connected particles 'offDiagonal' for each step;
                  connectedPart = obj.GetObjectConnectedParticles(oIdx,'offDiagonals');
                  obj.connectedStickyBeads = [obj.connectedStickyBeads;connectedPart];
                  end
                  %the particles connected can't be
                  %disconnected in the same step;
                  if ~isempty(row)
                    obj.connectedStickyBeads =   obj.connectedStickyBeads(1:(end-numel(stickyBeads(row))),:);
                  end
                  %add the probablity for each connected particles to test
                  %if they should be disconnected;
                  r = rand(size(obj.connectedStickyBeads,1),1);
                  disconnectInd = r<stickyBeadDisconnectProb;
                   for rIdx = 1:size(obj.connectedStickyBeads,1)
                          if disconnectInd(rIdx)
                              sprintf('%d and %d disconnected',obj.connectedStickyBeads(rIdx,1),obj.connectedStickyBeads(rIdx,2))
                              obj.DisconnectParticles(obj.connectedStickyBeads(rIdx,1),obj.connectedStickyBeads(rIdx,2))
                          end
                   end
                 obj.connectedStickyBeads = [];
        end
        
        function ConnectParticles(obj,particle1,particle2)% should move to ObjectInteractionManager
            % Connect two particles
            % particle1 and particle2 are number of particles in the general particle list   
            % if particle1 and particle2 are not members of the same object
            % the two containing objects are merged
            % Listeners to the event connectivityChange are informed only
            % if connectivity has changed upon merging of objects 
            
            % Check if the particles come from the same object 
            objNum = obj.map.GetObjectFromInd([particle1,particle2]);
%             obj2 = obj.map.GetObjectFromInd(particle2);
            
            % Change general connectivity % notifies listeners
            obj.connectivity(particle1,particle2) = true;
            obj.connectivity(particle2,particle1) = true;
              
            if objNum(1)~=objNum(2)
                obj.Merge((objNum))
                  disp('Merge')
            else
                disp('Same obj')
                notify(obj,'connectivityChange');
            end
            
        end
        
        function DisconnectParticles(obj,particle1,particle2)%should move to ObjectInteractionManager
            % Disconnect particle1 and particle2 
            % If the two particles are not members of the same objects, the
            % Split function is invoked. 
            
             % Check if the particles come from the same object 
             objNum = obj.map.GetObjectFromInd([particle1, particle2]);
%              obj2 = obj.map.GetObjectFromInd(particle2);             
             
             if objNum(1)~=objNum(2)
                 disp('particles are not connected')                               
             else
                 % find the members of the particles 
                 memb1 = obj.map.GetMemberFromInd(particle1);
                 memb2 = obj.map.GetMemberFromInd(particle2);
                 obj.connectivity(particle1, particle2) = false;
                 obj.connectivity(particle2, particle1) = false;
                  if memb1== memb2
                      % if it is the same member, update the connectivity map 
                      notify(obj, 'connectivityChange')                      
                      disp('disconnect beads of the same member')
                  else   
                    % Find connected  beads without nearest neighbors 
                    cParticles = obj.GetMembersConnectedParticles([memb1, memb2],'offDiagonals');
                    % get the connectivity between members
                    if isempty(cParticles)
                        % split the member from the object                        
                        membNumInObj = obj.map.GetMemberNumInObject(objNum(1),max(memb1,memb2));
                        obj.SplitMember(objNum(1),membNumInObj);
                        disp('split members ')
                    end                                        
                  end                 
             end
        end
                
        function Step(obj,objNum,dt)%TODO: find a uniform interface for composite and simple structures
            % Advance the objects one step in the simulation, apply forces
            % etc. 
            obj.particleDist = ForceManager.GetParticleDistance(obj.curPos);
            obj.particlesOnBoundary = [];
            cNb=0;             
            for oIdx = 1:numel(objNum)% for each object
                memberList = obj.map.GetObjectMembers(oIdx);% members of the objects
               
                if numel(memberList)==1
                        % get the indices for the members of the object 
                        beadInds     = obj.map.GetMemberInds(objNum(oIdx),1);
                        beadDistance = obj.particleDist(beadInds,beadInds);
                        % advance each member one step
                        obj.handles.chain(memberList(1)).Step(beadDistance,dt)   
                        % update the curPos list 
                        obj.curPos(beadInds,:)  = obj.handles.chain(memberList(1)).position.cur;

                else
                   % For composite object made of several sub-objects
                   connectivityMap     = obj.GetConnectivityMapAsOne(objNum(oIdx));
                   particleDistance    = obj.GetParticleDistance(objNum(oIdx)); 
                   [~,curMemberPos]    = obj.GetMembersPosition(objNum(oIdx));
                   springConst         = obj.GetSpringConstAsOne(objNum(oIdx));
                   minParticleDist     = obj.GetMinParticleDistAsOne(objNum(oIdx));
                   fixedParticleNum    = obj.GetFixedParticles(objNum(oIdx));
                   particlesBoundary   = obj.GetParticlesOnBoundary(objNum(oIdx));
                   
                   % split parameters to pass to the forceManager
                   par = obj.GetObjectParameters(objNum(oIdx));
                   par = par{1};
                   fp  = [par.forceParams];
                   
                   %TODO: work out fixed bead num for composite objects                  
                   newPos = ForceManager.ApplyCompositeInternalForces(curMemberPos,particleDistance,connectivityMap,...
                                                         [fp.springForce],[fp.bendingElasticityForce],...
                                                         springConst,[fp.bendingConst],...                                                         
                                                         minParticleDist,fixedParticleNum,particlesBoundary,dt);
                                                                                     
                % Deal the new pos to the object 
                   obj.DealCurrentPosition(oIdx,newPos);
                end
                %Deal with particles on the boundary for each step;
                obj.particlesOnBoundary = ...
                [obj.particlesOnBoundary, {obj.GetParticlesOnBoundary(objNum(oIdx)) + cNb}];
                %Update cumulative object indices 
                cNb = cNb+numel(obj.map.GetAllInds(oIdx)); 
            end 
          
            
            % Check for possible interaction between objects
%             obj.ObjectInteraction;
   
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
        
        function ObjectInteraction(obj)% move to objectInteractionManager
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
%                 o = 1:obj.map.count;% randperm(obj.numObjects);
%                 obj1 = obj.numObjects;
                
%                 rp = [1 100]; 
                if r>prob
                    
%                     pp = [1:32; 33:64]; 
                    rp = randperm(size(obj.connectivity,1));
%                     rp = [pp(1,rp(1)), pp(2,rp(2))];
%                      if obj.numObjects>1
%                           obj2 = obj.numObjects-1;
                    if  ~obj.connectivity(rp(1),rp(2))
                         obj.ConnectParticles(rp(1),rp(2));
                    % connect the head of the 2nd object and the tail of the
%                       1st one 
                      
% %                       get indices                       
%                       obj1Inds = (obj.map.GetAllInds(obj1));%obj.objectInds{obj.map.count};
%                       obj2Inds = (obj.map.GetAllInds(obj2));%obj.objectInds{obj.map.count-1};
%                       obj.connectivity(obj2Inds(1),obj1Inds(end)) = true;
%                       obj.connectivity(obj1Inds(end),obj2Inds(1)) = true;    
%                       % TODO: fix sorting in Merge function 
%                       obj.Merge(sort([obj1 obj2]));
%                       disp('merge')
                    end
%                      else
               
%                      end
                elseif r<(1-prob)
                    % get the list of connected particles, chose a pair and
                    % disconnect it 
                    [connected] = obj.GetObjectConnectedParticles(obj.map.count,'offDiagonals');
                    if ~isempty(connected)
                     % choose a pair 
%                      rn = randperm(size(connected,1));
%                      rp = connected(rn(1),:);
                     for cIdx = 1:size(connected,1)
                         obj.DisconnectParticles(connected(cIdx,1),connected(cIdx,2));
                     end
                    end
%                         m = obj.map.GetObjectCount(1:obj.numObjects);
%                          % split the biggest one 
%                         [numMem, oInd] = max(m);
%                         if numMem>1
%                         indsEnd   = obj.map.GetAllInds(oInd);%obj.objectInds{obj.map.count};
% %                         indsBend  = obj.map.GetMemberInds(oInd,2);%obj.objectInds{obj.map.count};
%                        obj.connectivity(indsEnd(1),indsEnd(end)) = false;
%                        obj.connectivity(indsEnd(end),indsEnd(1)) = false;%                       
%                        obj.SplitMember(oInd,obj.map.GetObjectCount(oInd))
%                          disp('split')
%                         end
              
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

                    
                   
              

            % ============================================
%             encounterDist = 0.1;
            % search for close beads and connect them, update the
            % connectivity accordingly 
%             obj.connectivity = obj.connectivity | obj.particleDist<encounterDist;
               end
        end
        
        function ConnectivityChangeListenerCallback(obj,varargin)
            % postSet - notify all registered listeners for connectivity
            % changes
            notify(obj,'connectivityChange')
        end
        
        function domainInd = GetDomainIndices(obj)
            % Give the domain number the particles belong to (all
            % particles)
            domainInd =  obj.domainInds;
        end
        
    end

methods (Static)

   function [row,col] = ListConnectivity(cMat)
            % find the connectivity of cMAt
           [row, col] = find(triu(cMat));
    end
end

end