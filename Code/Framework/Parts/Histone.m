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
        
        params = struct('numHistones',[],...% to be moved out
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
        
        function Step(obj,chainPosition,dt)
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
            fp                 = obj.params.forceParams;
            dirVec             = chainPos(obj.curPosVertex2,:)- chainPos(obj.curPosVertex1,:);
            diffusionForce     = ForceManager.GetDiffusionForce(fp.diffusionForce,obj.curPos,fp.diffusionConst,dt,[]);
            diffusionForce     = sqrt(sum(diffusionForce.^2,2));    
            % advence to a tentative location
            newTempPos         = obj.curPos+bsxfun(@times,diffusionForce,dirVec);% (obj.curPos(hIdx,:)-obj.prevPos(hIdx,:));

            obj.prevPos        = obj.curPos;
            obj.prevPosSlope   = obj.curPosSlope;
            for hIdx = 1:obj.params.numHistones
                         
                [obj.curPos(hIdx,:),vert1,vert2,obj.curPosSlope(hIdx)]= ...
                    obj.MoveOnChain(obj.curPos(hIdx,:),newTempPos(hIdx,:),obj.curPosVertex1(hIdx),obj.curPosVertex2(hIdx), chainPos,false);
                 if vert2 ~=obj.curPosVertex2 % particles changed edge                     
                  % update previous
                  obj.prevPosVertex1(hIdx) = obj.curPosVertex1(hIdx);
                  obj.prevPosVertex2(hIdx) = obj.curPosVertex2(hIdx); 
                  
                  % update current
                  obj.curPosVertex1(hIdx) = vert1; 
                  obj.curPosVertex2(hIdx) = vert2; 
                 end
            end
        end
        
        function UpdateHistonePositionOnChain(obj,chainPos)
            % Update the position of the histones after the chain has moved
%             obj.prevPos = obj.prevPos+ bsxfun(@times,obj.prevPosSlope,chainPos(obj.prevPosVertex2,:)- chainPos(obj.prevPosVertex1,:));
%                         obj.curPos = obj.curPos+ bsxfun(@times,obj.curPosSlope,chainPos(obj.curPosVertex2,:)- chainPos(obj.curPosVertex1,:));

             % updqte the current position on the new chain position
             curDirVec     = chainPos(obj.curPosVertex2,:)-chainPos(obj.curPosVertex1,:);
             obj.curPos = chainPos(obj.curPosVertex1,:)+bsxfun(@times,obj.curPosSlope,curDirVec);
             % update the previous chain position on the new chain position
             prevDirVec  = chainPos(obj.prevPosVertex2,:)-chainPos(obj.prevPosVertex1,:);
             obj.prevPos = chainPos(obj.prevPosVertex1,:)+bsxfun(@times,obj.prevPosSlope,prevDirVec);

%             for hIdx = 1:obj.params.numHistones
%                 dirVec = chainPos(obj.prevPosVertex2(hIdx,:),:)- chainPos(obj.prevPosVertex1(hIdx,:),:);
%                 % previous position
%                 obj.prevPos(hIdx,:) = chainPos(obj.prevPosVertex1(hIdx,:),:) +...
%                     obj.prevPosSlope(hIdx,:).*dirVec;
%                 % current position
%                 obj.curPos(hIdx,:)  = chainPos(obj.curPosVertex1(hIdx,:),:) + ...
%                     obj.curPosSlope(hIdx,:).*(chainPos(obj.curPosVertex2(hIdx,:),:)- chainPos(obj.curPosVertex1(hIdx,:),:));
% %                 % update the vertices 
% %                 obj.prevPosVertex1(hIdx) = obj.curPosVertex1(hIdx);
% %                 obj.prevPosVertex2(hIdx) = obj.curPosVertex2(hIdx);
% %                 [obj.curPosVertex1(hIdx), obj.curPosVertex2(hIdx)] = obj.FindPointOnPolygon(obj.curPos(hIdx,:),chainPos);% TODO: fix for ends 
%             end

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
                flag    = all(abs(t-t(1))<1e-10);
                vert1   = vertNum;
                vert2   = vertNum+1;
                vertNum = vertNum+1;
                if vertNum>numVert
                    flag = true;
                    vert1 = [];
                    vert2 = [];
                end
            end
        end
        
        function [newPos,vert1,vert2,posRatio] = MoveOnChain(prevPos,curPos,curPosVertex1, curPosVertex2,vertices,isCircular)%TODO: expand to Move on graph            
            % start from the prevPos between prevPosVertex1 and
            % prevPosVertex2 on the polygon with vertices in vertices
            % find the direction according to curPos-prevPos and relocate
            % the particle position to newPos between vert1 and vert2
            
            % project the curPos-prevPos on the segment of prevPos 
            pPath      = curPos-prevPos;
            % find the direction of motion (toward or away from polygon
            % start)            
            numPolyVert = size(vertices,1);% number of polygon vertices
            motionDir   = sum(pPath.*(vertices(curPosVertex2,:)-vertices(curPosVertex1,:)));
            motionDir   = sign(motionDir);
            pathLength  = sqrt(sum(pPath.^2));
            
            % start from prevPos and subtract the path length from polygon
            % length
            flag = false;
             vert1 = curPosVertex1;
             vert2 = curPosVertex2;
            if motionDir>0
               vert0 = vert2;
            else
               vert0 = vert1;
            end
             cumDist = sqrt(sum((prevPos-vertices(vert0,:)).^2));
             
            while ~flag
               
               pathReminder = pathLength-cumDist;
               flag         = pathReminder<0;
               if ~flag %  move to the next segment                   
                    vert1 = vert1+motionDir;
                    vert2 = vert2+motionDir;% TODO: insert reflection function                     
                    cumDist = cumDist+ sqrt(sum((vertices(vert1,:)-vertices(vert2,:)).^2));
               end
            end
                        
            edgeDir    = motionDir*(vertices(vert2,:)-vertices(vert1,:));  
            edgeLength = sqrt(sum((vertices(vert2,:)-vertices(vert1,:)).^2));
            posRatio   = (pathReminder+edgeLength)/edgeLength;% the position ratio between vertices
            
            if vert2==curPosVertex2 
                % if the particle did not switch edges                
                vert0 = vert1;               
            else
                % if the particle changed edges
             if motionDir>0
                 % motion in the 'positive' direction                
                 vert0 = vert1;                            
             else
                 % motion in the 'negative' direction              
                 vert0 = vert2;             
             end
            end
            % update the new position 
             newPos     = vertices(vert0,:)+ posRatio*edgeDir;
        end

    end
    
end