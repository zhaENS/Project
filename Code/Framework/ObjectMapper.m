classdef ObjectMapper<handle
    % Maps object numbers to their members indices in the general list of
    % ObjectManager    
    properties (SetObservable)
        count  = 0;
        object 
        members % sequence of all members
    end
    
    events
        countChange% change the number of objects 
    end
    
    methods
        
        function obj = ObjectMapper()
            % Class constructor
            
            % Add a callack to the countChange event
            addlistener(obj,'count','PostSet',@obj.CountChangeCallback);
        end
        
        function AddObject(obj,inds)
            % Add a member to the end of the list with the membes inds 
            obj.count          = obj.count+1;
            if ~isstruct(inds)
                objStruct                    = ObjectMapper.NewObjectStruct;
                objStruct.count              = 1;
                objStruct.members            = obj.count;% the index in the general list
                
                objStruct.inds.allInds       = [objStruct.inds.allInds,inds];
                objStruct.inds.memberInds{1} = inds;
                objStruct.inds.memberCount   = cellfun(@numel,objStruct.inds.memberInds);
                obj.members(obj.count)       = obj.count;  % permanent list of registered members
                obj.object                   = [obj.object;objStruct];

            else
                % Insert inds as a ready-made object structure
                % add the object structure 
                obj.object = [obj.object;inds];                
            end
        end
        
        function RemoveObject(obj,objList)
            % Remove an object from the list, do not alter the member list 
            
             stay = setdiff(1:obj.count,objList);             
             
             % Remove the objects
             obj.object = obj.object(stay);
                                       
             % Decrease object count
             obj.count  = obj.count-numel(objList);                          
             
        end
        
        function MergeObjects(obj,objList)
            % Merge two or more objects and their members into one
             % RemoveObjects->AddObject
%              objList = sort(objList);
             % copy data            
             memb     = ([obj.object(objList).members]);
             % make sure the list is sorted in ascending manner
             [memb,sIdx] = sort(memb,'ascend');             
             mCount   = numel(memb); % total number of members 
             oInds    = ([obj.object(objList).inds]);       % object inds           
             memCount = [oInds.memberCount];              % member count 
%              memCount = memCount(sInds);
             allInds  = ([oInds.allInds]);                  % indices from all members (FIFO)
             pMemb    = [oInds.memberInds];               % indices of each member
             
            
             % Remove all objects in objList
             obj.RemoveObject(objList);
             
             % Add the new object
             objStruct = ObjectMapper.NewObjectStruct;
             
             % fill the new struct
             objStruct.members          = memb;
             objStruct.count            = mCount;
             objStruct.inds.memberCount = memCount(sIdx);
             objStruct.inds.allInds     = sort(allInds);
             objStruct.inds.memberInds  = pMemb(sIdx);
             
             % add it to the end of the object list 
             obj.AddObject(objStruct);
                          
        end
        
        function BreakObject(obj,objNum)
            % Split an object to all its members 
            
            % sequential SplitMember
            while obj.object(objNum).count>1
                obj.SplitMember(objNum,obj.object(objNum).count)
            end
                        
        end
        
        function AddMember(obj,objNum,inds)
            % Add a new member to a specific object
            
            % AddObject -> MergeObjects   
            obj.AddObject(inds)
            obj.MergeObjects([objNum,obj.count])
        end
        
        function RemoveMember(obj,objNum,memberNum)
            % Remove a member from a specific object                      
            
            % remove the member's data from the object 
            obj.object(objNum).inds.memberInds  = obj.object(objNum).inds.memberInds(setdiff(1:obj.object(objNum).count,memberNum));
            obj.object(objNum).inds.allInds     = [obj.object(objNum).inds.memberInds{:}];
            obj.object(objNum).inds.memberCount = cellfun(@numel,obj.object(objNum).inds.memberInds);
            obj.object(objNum).count            = obj.object(objNum).count-1;
            obj.object(objNum).members          = setdiff(obj.object(objNum).members,obj.object(objNum).members(memberNum));            
            
        end
        
        function SplitMember(obj,objNum,memberNum)
            % Split a member of the object from the object 
             % RemoveMember-> AddObject     
            if obj.object(objNum).count>1
                             
             objStruct = obj.GetMemberAsObject(objNum,memberNum);
             obj.RemoveMember(objNum,memberNum);
             obj.AddObject(objStruct);
            end
        end
        
        function objStruct = GetMemberAsObject(obj,objNum,memberNum)
            % Wrap a member data as an object structure
            % the member number is relative to the order of appearance in the containing object  
            objStruct = ObjectMapper.NewObjectStruct;            
            % get member values 
            objStruct.count            = 1; % number of members
            objStruct.members          = obj.object(objNum).members(memberNum);          % member number            
            objStruct.inds.allInds     = obj.object(objNum).inds.memberInds{memberNum};  % all members indices
            objStruct.inds.memberInds  = obj.object(objNum).inds.memberInds(memberNum);  % member indices
            objStruct.inds.memberCount = obj.object(objNum).inds.memberCount(memberNum); % number of indices
        end
        
        function inds = GetAllInds(obj,objList)
            % Get all indices from different objects in one list
            objList = sort(objList);
            aInds   = [obj.object(objList).inds];
            inds    = ([aInds.allInds]);
        end
        
        function inds = GetMemberInds(obj,objNum,memberNum)
            % Get the indices of members of an object in one list
            % the member number is relative to the list appearing in the
            % containing object 
            memberNum = (memberNum);
            inds      = [obj.object(objNum).inds.memberInds{memberNum}];
        end
        
        function objNum = GetObjectFromInd(obj,ind)
            % get the object number as a fucntion of the index
            numInd = numel(ind);
            objNum = zeros(numInd,1);
            for iIdx = 1:numInd
             flag   = true;
%              next   = 1;
             objNum(iIdx) = 1;
             numObj = obj.count;
             while flag
                 % go over each object and test for inclusion of ind
                 inds = obj.GetAllInds(objNum(iIdx));
                 flag = ~ismember(ind(iIdx),inds);
                 if flag% if not there go to the next object 
                     objNum(iIdx) = objNum(iIdx)+1;
                     if objNum>numObj
                         error('the specified index is not in the list')
                     end
                 end
             end
            end
%             objNum = obj.allIndsToObj(ind);
        end
        
        function objList = GetObjectFromMember(obj,memberList)
            % get object number from member number (global member number)
            numMember = numel(memberList);
            objList   = zeros(numMember,1);
            for mIdx = 1:numMember
                flag  = false;
                objNum = 1; 
                while ~flag                     
                    flag = ismember(memberList(mIdx),obj.object(objNum).members);
                     if flag % if is a member of the object
                         objList(mIdx) = objNum; % record the object number 
                     else % if member of the object                         
                         objNum = objNum+1; % 
                         if objNum>obj.count % if not exceeds number of objects 
                             error('memberList is not on the list of members of any object');
                         end
                     end
                end                
            end
        end
        
        function memberNum = GetMemberFromInd(obj,indList)
            % get member number as a functio of the index
            memberNum = zeros(numel(indList),1);
            for iIdx = 1:numel(indList)
             objFlag = false;
             oIdx    = 1;             
             while ~objFlag 
                memFlag = false;
                mIdx    = 1;
                oCount  = obj.object(oIdx).count;
                 while ~memFlag 
                     membInds = obj.GetMemberInds(oIdx,mIdx);
                     memFlag  = ismember(indList(iIdx),membInds);
                     if memFlag 
                         memberNum(iIdx) = obj.object(oIdx).members(mIdx);
                         objFlag         = true;
                     else 
                         
                         mIdx = mIdx+1;% move to the next member 
                         
                         if mIdx>oCount% step out and move to the next object
                             memFlag = true;
                         end
                     end                                         
                 end                     
                     if oIdx==obj.count
                         objFlag = true;
%                          error('inds are not on the list')                     
                     else
                         oIdx = oIdx+1;
                     end 
              end
            end      
        end
        
        function count = GetObjectCount(obj, objNum)
            % Get the number of members for an object objNum
            count = [obj.object(objNum).count];
        end
        
        function count = GetMemberCount(obj,objNum,memberNum)
            % Number of indices for the member of an object 
            count = sum(obj.object(objNum).inds.memberCount(memberNum));
        end
        
        function members = GetObjectMembers(obj,objNum)
            % Get the member indices of object objNum
            members = obj.object(objNum).members;
        end
        
        function memberNum = GetMemberNumInObject(obj,objNum,memberList)
            % Get the member number according to the order of appearance in
            % object. The member number is as appears in the general member
            % list 
            [~,memberNum] = ismember(memberList,obj.object(objNum).members);
            
        end
        
        function CountChangeCallback(obj,varargin)
            % Notify listeners for the event of object count change 
            notify(obj,'countChange');
        end
    end
    
    methods(Static)
        
        function objectStruct = NewObjectStruct()
            % New structure for an object
            
            % objectStruct.count    - the number of members the object has 
            % objectStruct.members  - indices of members of object i (in
            % the fixed list)
            % objectStruct.inds.memberCount - number of indices of each member
            % objectStruct.inds.allInds     - all indices of all members
            % by order of appearance
            % objectStruct.inds.memberInds  - cell array of indices for 
            objectStruct = struct('count',0','members',[],...
                'inds', struct('memberCount',0,'allInds',[],'memberInds',cell(1)));
        end        
    end
    
end