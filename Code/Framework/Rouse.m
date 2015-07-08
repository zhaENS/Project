classdef Rouse<handle
    
    % TODO: add chain starting position options: linear, dense, loose, random,
    %       etc...
    % TODO: when beads are connected, make their stating position close.
    % make the starting position dense by restricting each bead to lay at a
    % distance from the center of mass (sequentially in creating the
    % initial beads)
    % TODO: add Reset option (so the class doesn't have to be initialized again)
    
    %========
    % Simulate the Rouse chain model
    % The system of differential equations
    % dR_n/dt = -dim*k*T/(s*b^2)*(2*R_n-R_{n+1} -
    % R_{n-1})+g_n(t)+BE_n(t)+LJ_n(t)
    % for the end beads
    % dR_1/dt = -dim*k*T/(s*b^2) *(R_1-R_2)+g_1(t) +BE_1(t)
    % dR_N/dt = -dim*k*T/(s*b^2) *(R_N-R_{N-1})+g_N(t)+BE_N(t)+LJ_N(t)
    
    % The function G_n(t) is the random fluctuation for the bead n, it is a
    % Gaussian with zero mean and 2*k*T/s std
    properties
        position
        connectionMap
        mobilityMatrices
        beadsDist
        relaxationTime %[sec] relaxation time of the first mode        
    end
    
    properties (SetObservable)
        params
        indsInParent
    end
    
    events
        step 
    end
    
    methods
        
        function obj = Rouse(rouseParams,indsInParent,parentHandle)
            % Class constructor
            % the chain parent object is the objectManager 
            % the inds in parent represent a list of integer in
            % objectManager;
            
            obj.params        = rouseParams;
            
            obj.indsInParent  = indsInParent; % the list of indices for position in ObjectManager
            if exist('parentHandle','var')
                % register a listener to the posChange event of
                % objectManager
                addlistener(parentHandle,'prevPosChange',@obj.UpdatePrevPosListenerCallback);
                addlistener(parentHandle,'curPosChange',@obj.UpdateCurPosListenerCallback);
                addlistener(parentHandle,'connectivityChange',@obj.UpdateConnectivityListenerCallback);
            end
            
            % Adjust input parameters to match the structure expected in
            % the class
            obj.SetInputParams;            
            
            obj.InitializeRouseStruct;
                        
        end        
        
        function UpdatePrevPosListenerCallback(obj,sourceObj,varargin)
            % pull the prevPosition from ObjectManager
            obj.position.prev(:,1:obj.params.dimension) = sourceObj.prevPos(obj.indsInParent,1:obj.params.dimension);
        end
        
        function UpdateCurPosListenerCallback(obj,sourceObj,varargin)
            % pull the curPosition from ObjectManager
            obj.position.cur(:,1:obj.params.dimension) = sourceObj.curPos(obj.indsInParent,1:obj.params.dimension);
        end
        
        function UpdateConnectivityListenerCallback(obj,sourceObj,varargin)
            % pull the connectivityMap from ObjectManager
            obj.connectionMap.map= sourceObj.connectivity(obj.indsInParent,obj.indsInParent);
            
        end
        
        function SetInputParams(obj)

            % expend the beta vector to match the number of beads 
            if numel(obj.params.beta)~=obj.params.numBeads && numel(obj.params.beta)~=1
                error('the number of beta values must match the number of beads')
            end             
            if numel(obj.params.beta)==1
                obj.params.beta = obj.params.beta*ones(obj.params.numBeads,1);
            end                        
        end
        
        function InitializeRouseStruct(obj)
            
            % Initialize positions
%             obj.position.cur             = randn(obj.params.numBeads,3);       
%             obj.position.cur(:,(1:3)>obj.params.dimension) = 0;
%             obj.position.prev            = randn(obj.params.numBeads,3);            
%             obj.position.prev(:,(1:3)>obj.params.dimension) = 0;

            obj.InitializeBeadConnectionMap;

        end        
        
        function CalculateRelaxationTime(obj)
            % the relaxation time of the first mode of the Rouse chain 
            d = sqrt(2*obj.params.diffusionConst*obj.params.dt);
            b = obj.params.b;
            N = obj.params.numBeads;
            p = 1;
            obj.relaxationTime = b^2/(12*d^2 *sin(p*pi/(2*N))^2);%[sec]
        end
        
        function InitializeBeadConnectionMap(obj,varargin)
            % Set the default linear connection 
            % define the Rouse matrix             
            cm = obj.params.connectedBeads;
            assert(all(cm(:)<=obj.params.numBeads),...
                'The monomers index to connect cannot exceed the number of monomers in the chain');
            
           % add the connection between distant beads
           cMap = logical(diag(ones(1,obj.params.numBeads-1),-1)+diag(ones(1,obj.params.numBeads-1),1));
           for bIdx = 1:size(cm,1)
               cMap(cm(bIdx,1),cm(bIdx,2))= true;
               cMap(cm(bIdx,2),cm(bIdx,1))= true;
           end
           
           obj.connectionMap.map                = logical(cMap);
                      
           [obj.connectionMap.indices.in.map]   = find(cMap); % linear indices of positions of monomer pairs connected. 
           [obj.connectionMap.indices.in.list(:,1), obj.connectionMap.indices.in.list(:,2)] = find(triu(cMap));
                       
        end        
         
        function UpdateLinearConnectivityMap(obj)
            % activated upon changing parameres of connectivityMap
          
             cMap                               = obj.connectionMap.map;
             [obj.connectionMap.indices.in.map] = find(cMap); % linear indices of positions of monomer pairs connected. 
             [inList(:,1), inList(:,2)]         = find(triu(cMap));
             obj.connectionMap.indices.in.list  = inList;
           
        end
        
        function [v,lambda] = GetRouseEig(obj)
            % Get the rouse matrix eigen values and eigen vectors
            v      = zeros(obj.params.numBeads);
            lambda = zeros(obj.params.numBeads,1);
            for bIdx = 1:obj.params.numBeads
                % the eigen vectors
                if bIdx ==1
                    v(:,bIdx)= 1/sqrt(obj.params.numBeads);
                else                    
                v(:,bIdx)    = sqrt(2/obj.params.numBeads)*cos(((1:obj.params.numBeads)-0.5)*pi*(bIdx-1)/(obj.params.numBeads))';
                end
                % the eigenvalues 
                lambda(bIdx) = 4*sin((bIdx-1)*pi/(2*(obj.params.numBeads))).^obj.params.beta(bIdx);% for the Rouse model set beta=2
                % in the present setting we allow the change of the
                % eigenvalues while the eigenvectors remain the same as in
                % the the original Rouse model 
            end
        end
        
        function SetInitialChainPosition(obj,domainHandler)%see comment
            % Start with a randomly stretched chain
            % the positions are initialized with zeros 
            % the secondInput domainHandler is the DomainHandler class. It
            % is an optional input, if inserted, the function uses it to
            % set the initial chain position such that all beads are inside
            % the domain, else it is randomly placed.
            
            %TODO: change to use the ForceManager
%             fp = domainHandler.params(obj.params.initializeInDomain).forceParams;
           
            % set initial position for the first beads
           flag = false;
           dimInds = (1:3)>obj.params.dimension;% dimensions not included
           exDimDomVals = domainHandler.params(obj.params.initializeInDomain).domainCenter(dimInds); 
           while ~flag
            obj.position.prev(1,:)       = 0.2*obj.params.b*randn(1,3);
            obj.position.prev(1,dimInds) = exDimDomVals;
            flag = domainHandler.InDomain(obj.position.prev(1,:),obj.params.initializeInDomain);      
           end
           
            if exist('domainHandler','var')
              if isempty(obj.params.beadsOnBoundary)
                fp = domainHandler.params(obj.params.initializeInDomain).forceParams;
                % The bead positions                
                for bIdx = 2:obj.params.numBeads         
                    inDomain = false;
                    while ~inDomain  
                     dx          = sqrt(2*fp.diffusionConst)*randn(1,3);
%                      dx          = sqrt(2*fp.diffusionConst)*randn(1,3);
                     dx(dimInds) = exDimDomVals;
                     tempPos     = obj.position.prev(bIdx-1,:)+dx;
                     inDomain    = domainHandler.InDomain(tempPos,obj.params.initializeInDomain);     
                    end                                          
                    obj.position.prev(bIdx,:)= tempPos;
                    
                end
                obj.position.prev(obj.params.fixedBeadNum,:) = obj.params.fixedBeadsPosition;
                                
                
                else % if there are beads constrained to lay on the boundary 
                    
                    % 1. diffuse on the boundary from an initial point to
                    % set position for the constrained beads
                    % 2. iteratively pass a Brownian bridge between these
                    % points
                    numSteps   = range(obj.params.beadsOnBoundary); % this serves as the number of steps to walk on the boundary 
                    % get initial points 
                    initialPoint = domainHandler.GetRandomBoundarySample(1);
                    % diffuse from the initial point numSteps
                    % draw two angles from a normal wrapped distribution and advance accordingly
                    
             end
        
            
            else                
                % The bead positions
               r = randn(obj.params.numBeads-1,obj.params.dimension);
               obj.position.prev = [obj.position.prev(1,:); cumsum(r)];
            end
            
            obj.position.cur = obj.position.prev;
        end                        
                                        
        function Step(obj,beadDist,dt)
            
            % Apply forces on the beads to get the new bead position    
            forceParams = obj.params.forceParams;
            newPos      = ForceManager.ApplyInternalForces(obj.position.cur,beadDist,obj.connectionMap.map,...  
                                              forceParams.springForce,forceParams.bendingElasticityForce,...
                                              forceParams.springConst,forceParams.bendingConst,...   
                                              forceParams.bendingAffectedParticles,forceParams.bendingOpeningAngle,...
                                              forceParams.minParticleEqDistance,obj.params.fixedBeadNum,dt);

            obj.position.cur(:,1:obj.params.dimension) = newPos(:,1:obj.params.dimension);  % update current position % the positions are updated after the domain has exerted its force on the chain                      
        end
        
        function SetPrevBeadPosition(obj,pos)% obsolete, externaly used
             % Save previous position as the currentPosition
             if ~exist('pos','var')
               obj.position.prev = obj.position.cur; 
             else
               obj.position.prev = pos;
             end                          
        end                                   
        
        function ConnectBeads(obj,bead1,bead2)
            % Check if the pair is already in the list 
            flag = false;
            if isempty(obj.params.connectedBeads)
                flag = true;
            else
                % check if the pair is already on the list 
               b = any(all(bsxfun(@eq,[bead1,bead2],[obj.params.connectedBeads; fliplr(obj.params.connectedBeads)]),2));
            if ~b
              flag = true;
            end             
            end
            
            if flag 
                % Add the pair entry to the end of the connectivity list 
                obj.params.connectedBeads(end+1,:) = [bead1 bead2];
                obj.connectionMap.map(bead1,bead2) = true;
                obj.connectionMap.map(bead2,bead1) = true;
                
                % Update linear connectivity list
                obj.UpdateLinearConnectivityMap
            end
        end
        
        function DisconnectBeads(obj,bead1,bead2)
            
            % zero out the entries in the connectivity matrix 
            obj.connectionMap.map(bead1,bead2) = false;
            obj.connectionMap.map(bead2,bead1) = false;
             
            if ~isempty(obj.params.connectedBeads)
                
            % Remove the pair entry in the list 
            b1 = all(bsxfun(@eq,[bead1,bead2],[obj.params.connectedBeads]),2);
            b2 = all(bsxfun(@eq,[bead2,bead1],[obj.params.connectedBeads]),2);
            obj.params.connectedBeads((b1|b2),:) =[]; 
            
            % Update linear indices
            obj.UpdateLinearConnectivityMap
            end
        end
    end
        
end

