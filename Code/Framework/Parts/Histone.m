classdef Histone<handle
    % a class for a Histone object
    % histones are attached to chains by default and move on them by 1D
    % diffusion, the position in space is determined by the position of the
    % chain
    % the new chain position is passed to histone, the histone is then
    % projected onto the new position between the vertices it was on in the
    % previous step and then it is moved to its new position
    
    % TODO: fix FindPointOnPolygon to allow reflection at the ends 
    % TODO: omit the main loop to work on all the coordinates as a matrix 
    % TODO: expand MoveOnPolygon to MoveOnGraph
    
    properties
        curPos
        curPosVertex1
        curPosVertex2
        curPosSlope % position between vertices- used for transformation
        prevPos
        prevPosSlope % position between vertices- used for transformation
        prevPosVertex1
        prevPosVertex2
        
        params = struct('numHistones',[],...
                        'dt',[],...
                        'minHistoneDist',[],...
                        'diffusionForce',[],...
                        'diffusionConst',[],...                        
                        'ljForce',[],...
                        'LJPotentialWidth',[],...
                        'LJPotentialDepth',[],...
                        'morseForce',[],...
                        'morsePotentialDepth',[],...
                        'morsePotentialWidth',[],...
                        'morseForceType',[],...
                        'forceParams',ForceManagerParams);
                        
    end
    
    methods
        
        function obj = Histone(varargin)
            % class constructor 
            obj.ParseInputParams(varargin);
                        
        end
        
        function Initialize(obj,initialChainPosition)
            % initialize histone position
            
            % choose a random starting point on the chain
            % choose a bond
            chainPosition = initialChainPosition;
            % choose a random bond
            for hIdx = 1:obj.params.numHistones
                bondNum  = randperm(size(chainPosition,1)-1);
                
                s1 = chainPosition(bondNum(1),:);
                s2 = chainPosition(bondNum(1)+1,:);
                r  = sort(rand(1,2));
                % linearily interpolate the bead position and choose a random
                % point on the bond
                
                obj.prevPos(hIdx,:)        = s1+r(1)*(s2-s1);
                obj.prevPosVertex1(hIdx,:) = bondNum(1);
                obj.prevPosVertex2(hIdx,:) = bondNum(1)+1;
                obj.prevPosSlope(hIdx,:)   = r(1);
                
                obj.curPos(hIdx,:)        = s1+r(2)*(s2-s1);
                obj.curPosVertex1(hIdx,:) = bondNum(1);
                obj.curPosVertex2(hIdx,:) = bondNum(1)+1;
                obj.curPosSlope(hIdx,:)   = r(2);
            end
        end
        
        function Step(obj,chainPosition)
            % move by 1D diffusion
            % project onto new chain position
            chainPos          = chainPosition;
%             particleDistances = ForceManager.GetParticleDistance(obj.curPos);
            
            % translate the prev and cur pos onto the new chain positon
            % (after it moved)
%                         obj.prevPos = obj.prevPos+ bsxfun(@times,obj.prevPosSlope,chainPos(obj.prevPosVertex2,:)- chainPos(obj.prevPosVertex1,:));
%                         obj.curPos = obj.curPos+ bsxfun(@times,obj.curPosSlope,chainPos(obj.curPosVertex2,:)- chainPos(obj.curPosVertex1,:));
            %             % Transform curPos to distance on 1D line
            %
            %
            %             [newTempPos, diffusionForce,ljForce,~] = ForceManager.ApplyExternalForces(obj.curPos,...
            %                                                           particleDistances,obj.params.diffusionConst,...
            %                                                           ljForce,diffusionForce,morseForce,...
            %                                                           LJPotentialWidth,LJPotentialDepth,...
            %                                                           morsePotentialDepth, morsePotentialWidth,morseForceType,...
            %                                                           minParticleDist,fixedParticleNum,dt);
            
            % Update the histone position on the new chain position 
            obj.UpdateHistonePositionOnChain(chainPos);
            
            for hIdx = 1:obj.params.numHistones
              
                 dirVec = chainPos(obj.prevPosVertex2(hIdx,:),:)- chainPos(obj.prevPosVertex1(hIdx,:),:);
%                 % previous position
%                 obj.prevPos(hIdx,:) = chainPos(obj.prevPosVertex1(hIdx,:),:) +...
%                     obj.prevPosSlope(hIdx,:).*dirVec;
%                 % current position
%                 obj.curPos(hIdx,:)  = chainPos(obj.curPosVertex1(hIdx,:),:) + ...
%                     obj.curPosSlope(hIdx,:).*(chainPos(obj.curPosVertex2(hIdx,:),:)- chainPos(obj.curPosVertex1(hIdx,:),:));
                
                % get diffusion force
                diffusionForce = ForceManager.GetDiffusionForce(obj.params.diffusionForce,obj.curPos(hIdx,:),obj.params.diffusionConst,obj.params.dt,[]);
                
                % project the diffusion force on the current chain segment
                % to get the position                 
               
                diffusionForce      = sqrt(sum(diffusionForce.^2)); 
                newTempPos          = obj.curPos(hIdx,:)+diffusionForce.*dirVec;% (obj.curPos(hIdx,:)-obj.prevPos(hIdx,:));
                newPos              = obj.MoveOnChain(obj.curPos(hIdx,:),newTempPos,chainPos,false);
                obj.prevPos(hIdx,:) = obj.curPos(hIdx,:);
                obj.curPos(hIdx,:)  = newPos;
                % update the vertices 
                obj.prevPosVertex1(hIdx) = obj.curPosVertex1(hIdx);
                obj.prevPosVertex2(hIdx) = obj.curPosVertex2(hIdx);
                [obj.curPosVertex1(hIdx), obj.curPosVertex2(hIdx)] = obj.FindPointOnPolygon(obj.curPos(hIdx,:),chainPos);% TODO: fix for ends 
                
            end
        end
        
        function UpdateHistonePositionOnChain(obj,chainPos)
            % Update the position of the histones after the chain has moved
%             obj.prevPos = obj.prevPos+ bsxfun(@times,obj.prevPosSlope,chainPos(obj.prevPosVertex2,:)- chainPos(obj.prevPosVertex1,:));
%                         obj.curPos = obj.curPos+ bsxfun(@times,obj.curPosSlope,chainPos(obj.curPosVertex2,:)- chainPos(obj.curPosVertex1,:));
            for hIdx = 1:obj.params.numHistones
                dirVec = chainPos(obj.prevPosVertex2(hIdx,:),:)- chainPos(obj.prevPosVertex1(hIdx,:),:);
                % previous position
                obj.prevPos(hIdx,:) = chainPos(obj.prevPosVertex1(hIdx,:),:) +...
                    obj.prevPosSlope(hIdx,:).*dirVec;
                % current position
                obj.curPos(hIdx,:)  = chainPos(obj.curPosVertex1(hIdx,:),:) + ...
                    obj.curPosSlope(hIdx,:).*(chainPos(obj.curPosVertex2(hIdx,:),:)- chainPos(obj.curPosVertex1(hIdx,:),:));
            end


        end
        
        function [newPos,vert1,vert2] = MoveOnChain(obj,prevPos,curPos,vertices,isCircular)%TODO: expand to Move on graph
            %             A = [vertices(1,1:2),1; vertices(2,1:2),1 ; prevPos(1,1:2),1];
            %             d = det(A);
            %             if d>eps
            %                 error('prevPos is not on the first edge')
            %             end
            % Get the exct position and trim the polygon such that the polygon and the
            % path start from the same positon
            
            % Length of the path
            l       = sqrt(sum((curPos-prevPos).^2));
            flag    = true;
            
            % find the location of prevPos on the chain
            
            %             chainLength   = cumsum(sqrt(sum(vertices.^2,2)),1);
            [vert1,vert2] = obj.FindPointOnPolygon(prevPos,vertices);
            if isempty(vert1) || isempty(vert2)
                error('particle not on the polygon')
            end
            %             vert1      = 1;
            %             vert2      = 2;
            polyLength = 0;
            numVert    = size(vertices,1);
            while flag
                % iteratively calculate the length of the polygon
                % and subtract it from the length of the path
                r          = sqrt(sum((vertices(vert2,:)-vertices(vert1,:)).^2));  % length of segment
                polyLength = polyLength + r;% cumulative length of polygon
                d          = l-polyLength; % difference between path length and cumulative polygon length
                if d<0 % if the difference is negative- stop
                    flag = false;
                else
                    vert1 = vert2;
                    
                    if vert2==numVert
                        if isCircular
                            vert2 = 1;
                        else
                            %  quit with no result
                            error('the point is not on the specified polygon')
                        end
                    else
                        vert2 = vert2+1;% advance the index one more
                    end
                end
            end
            
            % the actual point
            a      = vertices(vert1,:);
            b      = vertices(vert2,:);
            dirVec = (b-a)./sqrt(sum((b-a).^2));
            
            t      = (r+d)/r;
            % get the new point
            newPos = a+t*(dirVec);
            
            polyVert = vertices;
            if isCircular
                polyVert = [polyVert; vertices(1,:)];
            end
        end
        
        function ParseInputParams(obj,varargin)
            % parse the input parameters
            p = varargin{:};
            if mod(numel(p),2)~=0
                error('name-value pair must be inserted')
            end
            
            for pIdx = 1:numel(p)/2
                obj.params.(p{2*pIdx-1}) = p{2*pIdx};
            end
            
        end
    end
    
    methods (Static)
        
        function [vert1,vert2]= FindPointOnPolygon(point,vertices)
            % find a point on a polygon between vert1 and vert2
            flag    = false;
            vertNum = 1;
            numVert = size(vertices,1);
            if numVert <2
                error('polygon must contain at least 2 vertices')
            end
            
            while~flag
                t       = (point - vertices(vertNum,:))./(vertices(vertNum+1,:)-vertices(vertNum,:));
                flag    = all(abs(t-t(1))<1e-9);
                vert1   = vertNum;
                vert2   = vertNum+1;
                vertNum = vertNum+1;
                if vertNum>=numVert
                    flag = true;
                    vert1 = [];
                    vert2 = [];
                end
            end
        end
        
    end
    
end