classdef ObjectMapper<handle
    % maps object numbers to their members indices in the general list of
    % ObjectManager    
    properties
        count  = 0;
        object 
        members % sequence of all members
    end
    
    methods
        function obj = ObjectMapper()
            % Class constructor
        end
        
        function AddObject(obj,inds)
            % Add a member to the end of the list
            obj.count          = obj.count+1;
            % obj.object(i).count    - the number of members the object has 
            % obj.object(i).members  - indices of members of object i (in
            % the fixed list)
            % obj.object(i).inds.memberCount - number of indices of each member
            % obj.object(i).inds.allInds     - all indices of all members
            % by order of appearance
            % obj.object(i).inds.memberInds  - cell array of indices for 
            
            objStruct                  = ObjectMapper.NewObjectStruct;                        
            objStruct.count            = 1;
            objStruct.members          = obj.count;% the index in the general list
            objStruct.inds.memberCount = objStruct.inds.memberCount+1;
            objStruct.inds.allInds     = [objStruct.inds.allInds,inds];
            objStruct.inds.memberInds{objStruct.inds.memberCount}= inds;
             
            
            obj.members(obj.count) = obj.count;
          
             
            obj.object             = [obj.object;objStruct];
        end
        
        function RemoveObject(obj,objList)
            % Remove an object from the list, do not alter the member list 
            
             stay = setdiff(1:obj.count,objList);
             % Get the list of members from objects in obj list 
%              memToRemove = [obj.object(objList).members];
%              obj.members = setdiff(obj.members,memToRemove);
             
             % Remove the objects
             obj.object = obj.object(stay);
             
             % Decrease object count
             obj.count  = obj.count-numel(objList);
        end
        
        function MergeObjects(obj,objList)
            % Merge two or more objects and their members into one
             
             memb     = [obj.object(objList).members];
             mCount   = sum([obj.object(objList).count]); 
             oInds    = [obj.object(objList).inds]; % object inds
             memCount = sum([oInds.memberCount]); % member count 
             allInds  = [oInds.allInds];
             pMemb    = [oInds.memberInds];
             
             % remove all members 
             obj.RemoveObject(objList);
             
             objStruct = ObjectMapper.NewObjectStruct;
             
             % fill the new struct
             objStruct.members = memb;
             objStruct.count   = mCount;
             objStruct.inds.memberCount = memCount;
             objStruct.inds.allInds     = allInds;
             objStruct.inds.memberInds  = pMemb;
             
             % add it to the end of the object list 
             obj.count = obj.count+1;
             obj.object = [obj.object;objStruct];
        end
        
        function BreakObject(obj,objNum)
            % split the object to all its members 
            
            % sequential SplitMember
            for mIdx=1:obj.object(objNum).count
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
            % remove a member from a specific object
            
            % SplitMember -> RemoveObject             
            obj.SplitMember(objNum,memberNum);
            obj.RemoveObject(obj.count);
        end
        
        function SplitMember(obj,objNum,memberNum)
            % Split a member of the object from the object 
            
            if obj.object(objNum).count>1
            % RemoveMember-> AddObject
            objStruct = ObjectMapper.NewObjectStruct;
            
            % get member values 
            objStruct.count            = 1;
            objStruct.members          = obj.object(objNum).members(memberNum);
            objStruct.inds.memberCount = 1;
            objStruct.inds.allInds     = obj.object(objNum).inds.memberInds{memberNum};
            objStruct.inds.memberInds  = obj.object(objNum).inds.memberInds(memberNum);
            
            % modify the object
            obj.object(objNum).count = obj.object(objNum).count-1;
            obj.object(objNum).members = setdiff(obj.object(objNum).count,objStruct.members);
            
            % Add the new object to the end of the list 
             obj.object = [obj.object;objStruct];
            end
        end
        
        function inds = GetAllInds(obj,objList)
            inds = [obj.object(objList).inds.allInds];
        end
        
        function inds = GetMemberInds(obj,objNum,memberNum)
            inds = obj.object(objNum).inds.memberInds{memberNum};
        end
        
        function count = GetObjectCount(obj, objNum)
            % Get the total number of indices of an object
            count = obj.GetMemberCount(objNum,1:obj.object(objNum).count);
        end
        
        function count = GetMemberCount(obj,objNum,memberNum)
            % Count of the member of an object 
            count = numel([obj.object(objNum).inds.memberInds{memberNum}]);
        end
        
        function members = GetObjectMembers(obj,objNum)
            members = obj.object(objNum).members;
        end
    end
    
    methods(Static)
        
        function objectStruct = NewObjectStruct()
            % new structure for an object
            objectStruct = struct('count',0','members',[],...
                'inds', struct('memberCount',0,'allInds',[],'memberInds',cell(1)));
        end        
    end
    
end