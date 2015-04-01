classdef Rouse<handle
    
    % TODO: add chain starting position options: linear, dense, loose, random,
    %       etc...
    % TODO: when beads are connected, make their stating position close.
    % make the starting position dense by restricting each bead to lay at a
    % distance from the center of mass (sequentially in creating the
    % initial beads)
    % TODO: add Reset option (so the class doesn't have to be initialized again)
    % TODO: add the possibility to change the Raouse matrix's eigenvalues,
    %       making the model a beta model
    %
    
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
        params
        handles
        positions                       
        step 
        time
        simulationTime
        distributions% should e removed  
        connectionMap

        runSimulation
        mobilityMatrices
        
        relaxationTime %[sec] relaxation time of the first mode
        
        beadsDist % move to ForceManager
        forces    % should move to ForceManager
        dimNames = {'x','y','z'}; % should be discarded, move to matrix representation 
    end
    
    properties (Access=private)
        noise
        noiseIndex = 0;
        numNoisePoints
    end
    
    events
        getForce% next step event 
        
    end
    
    methods
        function obj = Rouse(rouseParams)
            % class constructor
            obj.runSimulation = true;
            obj.params        = rouseParams;
            
            % Initialize the forces acting on the chain 
            obj.handles.classes.forceManager = ForceManager();  
            obj.handles.classes.forceManager.springForce       = obj.params.springForce;
            obj.handles.classes.forceManager.lennardJonesForce = obj.params.lennardJonesForce;
            obj.handles.classes.forceManager.diffusionForce    = obj.params.diffusionForce;
            obj.handles.classes.forceManager.bendingElasticityForce = obj.params.bendingElasticityForce;
            % Adjust input parameters to match the structure expected in
            % the class
            obj.SetInputParams
                        
            obj.InitializeRouseStruct
            
            
        end        
        
        function SetInputParams(obj)% clean up 
%             if mod(numel(paramsIn),2)~=0
%                 error('value pair input must come in pairs')
%             end
%             for pIdx = 1:2:numel(paramsIn);
%                 obj.params.(paramsIn{pIdx})= paramsIn{pIdx+1};
%             end
            % expend the beta vector to match the number of beads 
            if numel(obj.params.beta)~=obj.params.numBeads && numel(obj.params.beta)~=1
                error('the number of beta values must match the number of beads')
            end             
            if numel(obj.params.beta)==1
                obj.params.beta = obj.params.beta*ones(obj.params.numBeads,1);
            end
                
%             % Update noise STD
%             obj.params.noiseStd          = sqrt(2*obj.params.diffusionConst*obj.params.dt);
%             % update spring const.
%             obj.params.springConst       = -obj.params.dimension*obj.params.diffusionConst*obj.params.dt/(obj.params.b^2);
        end
        
        function InitializeRouseStruct(obj)
            
            obj.time(1)                         = 0;
            obj.step                            = 1;
            obj.simulationTime.total            = 0; % in seconds 
            obj.simulationTime.meanStepTime     = 0; % in seconds
            
            obj.positions.beads.cur.x           = zeros(obj.params.numBeads,1);
            obj.positions.beads.cur.y           = zeros(obj.params.numBeads,1);
            obj.positions.beads.cur.z           = zeros(obj.params.numBeads,1);
            
            obj.positions.beads.prev.x          = zeros(obj.params.numBeads,1);
            obj.positions.beads.prev.y          = zeros(obj.params.numBeads,1);
            obj.positions.beads.prev.z          = zeros(obj.params.numBeads,1);            
            
            obj.positions.springs.x             = sparse(1:obj.params.numBeads,1:obj.params.numBeads,0); 
            obj.positions.springs.y             = sparse(1:obj.params.numBeads,1:obj.params.numBeads,0);  
            obj.positions.springs.z             = sparse(1:obj.params.numBeads,1:obj.params.numBeads,0); 
            obj.positions.springs.lengths       = sparse(1:obj.params.numBeads,1:obj.params.numBeads,0);
            
            % The angle between springs is defined as a 3D array
            % but is represented as N 2D arrays of sparse matrices 
            obj.positions.springs.angleBetweenSprings = sparse(repmat(1:obj.params.numBeads,1,obj.params.numBeads),1:obj.params.numBeads^2,pi);
            obj.beadsDist                       = sparse(1:obj.params.numBeads,1:obj.params.numBeads,0);
            obj.connectionMap.beadTriplets      = [];
            
            % Initialize forces struct
            obj.forces.springs      = zeros(obj.params.numBeads);
            obj.forces.lennardJones = struct('x',zeros(obj.params.numBeads,1),...
                                             'y',zeros(obj.params.numBeads,1),...
                                             'z',zeros(obj.params.numBeads,1));
            obj.forces.bending = struct('x',zeros(obj.params.numBeads,1),...
                                        'y',zeros(obj.params.numBeads,1),...
                                        'z',zeros(obj.params.numBeads,1));        
                                
            % Initialize structures 
%             obj.SetBeadsMobilityMatrices;  
            obj.SetBeadConnectionMap
            obj.SetInitialChainPosition; % should be moved out of the initialization process

        end
                
        function CalculateRelaxationTime(obj)
            % the relaxation time of the first mode of the Rouse chain 
            d = sqrt(2*obj.params.diffusionConst*obj.params.dt);
            b = obj.params.b;
            N = obj.params.numBeads;
            p = 1;
            obj.relaxationTime = b^2/(12*d^2 *sin(p*pi/(2*N))^2);%[sec]
        end
        
        function SetBeadConnectionMap(obj)
            % Set the default linear connection 
            % define the Rouse matrix 
            obj.mobilityMatrices.rouse = RouseMatrix(obj.params.numBeads);                        
            
%             [v,lambda]           = obj.GetRouseEig;
            % if the lambdas are changed, then the connection map is
            % changed. 

            cm = obj.params.connectedBeads;
            assert(all(cm(:)<=obj.params.numBeads),'The monomers index to connect cannot exceed the number of monomers in the chain')            
            % recalculate the connection map 
%             connectionMatrix = v*diag(lambda)*v';
%             connectionMatrix = connectionMatrix.*(abs(connectionMatrix)>1e-12);
%             for bIdx = 1:size(cm,1)
%                 connectionMatrix(cm(bIdx,1),cm(bIdx,2))  = -1;
%                 connectionMatrix(cm(bIdx,2),cm(bIdx,1))  = -1;
%             end
            
%             obj.connectionMap.connectionMatrix = connectionMatrix; % used for multiplication of beads location 
%             cMap = logical((connectionMatrix-diag(diag(connectionMatrix)))~=0);
%             cMap(~cMap) = NaN;
           % add the connection between distant beads
           cMap = logical(diag(ones(1,obj.params.numBeads-1),-1)+diag(ones(1,obj.params.numBeads-1),1));
           for bIdx = 1:size(cm,1)
               cMap(cm(bIdx,1),cm(bIdx,2))= true;
               cMap(cm(bIdx,2),cm(bIdx,1))= true;
           end
           
           obj.connectionMap.map                = cMap;
           [obj.connectionMap.indices.in.map]   = find(cMap); % linear indices of positions of monomer pairs connected. 
           [obj.connectionMap.indices.in.list(:,1), obj.connectionMap.indices.in.list(:,2)] = find(triu(cMap));
%            [obj.connectionMap.indices.out.map]  =  ~cMap;% linear indices of positions of monomer pairs connected. 
%            [obj.connectionMap.indices.out.list(:,1), obj.connectionMap.indices.out.list(:,2)]  = find(obj.connectionMap.indices.out.map); 
%            cMapNz = cMap;
%            cMapNz(isnan(cMapNz))=0;
%            obj.connectionMap.numConnections = sum(cMapNz,2);% number of connection each bead has 
% 
%            % Set angle between springs as a triplet of beads. The middle
%            % bead represents the middle node in the triplet defining the
%            % angle. The convention is row number, column number, height
%            % number. To represent an angle between bead i,j,k. start from
%            % row i to column j and go to the depth layer k
%            
%            sMap = zeros(obj.params.numBeads,obj.params.numBeads,obj.params.numBeads);
%            % first register the linear connections 
%            for bIdx = 2:obj.params.numBeads-1
%                sMap(bIdx-1,bIdx,bIdx+1) = 1;
%                obj.connectionMap.beadTriplets(end+1,:) = [bIdx-1,bIdx,bIdx+1];
%            end
%            
%            % next, register the non-trivial connections
% %            cMap1 = triu(cMap);
% %            cMap1(isnan(cMap1))=0;
% %            cMap1(isnan(cMap1))=0;
%            for bIdx = 1:obj.params.numBeads
%                % there are either 2 or more non-zero elements in each row.
%                allPairs = [];
%                 f= find(cMap(:,bIdx)==1);
%                 for fIdx = 1:numel(f);
%                     g = find(cMap(f(fIdx),:)==1);
%                     g= setdiff(g,bIdx);
%                     for gIdx = 1:numel(g)
%                         allPairs(end+1,:) = [min([bIdx g(gIdx)]),f(fIdx),max([bIdx,g(gIdx)])];
%                     end
%                 end
% 
% %                 if numel(f)>2
% %                     allPairs = [f,bIdx*ones(numel(f),1),circshift(f,1)];
% %                 elseif numel(f)==2
% %                     allPairs = [f(1),bIdx,f(2)];                    
% %                 else 
% %                     allPairs = [];
% %                 end
%                 % sort the indices such that the smaller index is on the
%                 % left
% %                 for pIdx = 1:size(allPairs,1)
% %                     if allPairs(pIdx,1)>allPairs(pIdx,3)
% %                         allPairs(pIdx,:) = fliplr(allPairs(pIdx,:));
% %                     end                        
% %                 end
% %                 allPairs = unique(allPairs,'rows');
%                 for pIdx = 1:size(allPairs,1)
%                     % make a symmetric matrix in 3D
%                     sMap(allPairs(pIdx,1),allPairs(pIdx,2),allPairs(pIdx,3))= 1;
%                     sMap(allPairs(pIdx,3),allPairs(pIdx,2), allPairs(pIdx,1))=1;% check                     
%                     obj.connectionMap.beadTriplets(end+1,:) = allPairs(pIdx,:);
%                 end
%               obj.connectionMap.angles = sMap;  
%            end
%             obj.connectionMap.beadTriplets= unique(obj.connectionMap.beadTriplets,'rows');
%            % define the linear indices in the array of bead triplets to
%            % save running time later in the calculation of angle between
%            % beads. angles is the angle at the middle bead (bead2) so we
%            % have bead1ToBead2 and bead3ToBead2 to define the edges of the
%            % angle
%            bt     = obj.connectionMap.beadTriplets;
%            for btIdx = 1:numel(bt(:,1))
%             obj.connectionMap.indices.in.linear.bead1ToBead2(btIdx,1)  = sub2ind([obj.params.numBeads,obj.params.numBeads],min([bt(btIdx,1),bt(btIdx,2)]),max([bt(btIdx,1),bt(btIdx,2)]));
%             obj.connectionMap.indices.in.linear.bead3ToBead2(btIdx,1)  = sub2ind([obj.params.numBeads,obj.params.numBeads],min([bt(btIdx,2),bt(btIdx,3)]),max([bt(btIdx,2),bt(btIdx,3)]));
%            end
%            % bead triplets represents the indices in the N 2D sparse matrix
%            % of bead connections, row for bead1, bead2 is represented N
%            % times in N*2D matrix, bead3 is the column Nth block meaning
%            % size(matrix,2)/numBeads.                       
%            obj.connectionMap.indices.in.linear.beadTriplets  = (obj.connectionMap.beadTriplets(:,3)-1)*obj.params.numBeads^2+obj.connectionMap.indices.in.linear.bead1ToBead2;
%            k = sub2ind([obj.params.numBeads, obj.params.numBeads],obj.connectionMap.indices.in.list(:,1),obj.connectionMap.indices.in.list(:,2));
%            obj.connectionMap.indices.in.linear.beadPairs  = k;% pair of beads connected (linear indices of the matrix size [numBead X numBeads])
%            [c,~,ic] = unique(obj.connectionMap.beadTriplets(:,2),'R2012a');                
%            obj.connectionMap.indices.uniqueBeadTriplets.uniqueList = c;
%            obj.connectionMap.indices.uniqueBeadTriplets.uniqueListIndices = ic;
%            n       = 1:numel(obj.connectionMap.beadTriplets(:,2));
%            k       = [c(ic),n'];  
%            obj.connectionMap.indices.uniqueBeadTriplets.uniqueIndicesMap = k;
%            n = [1:obj.params.numBeads,1:abs(numel(obj.connectionMap.beadTriplets(:,2))-obj.params.numBeads)];
%            k       = obj.connectionMap.indices.uniqueBeadTriplets.uniqueIndicesMap;
%            obj.connectionMap.indices.uniqueBeadTripletsLinear = sub2ind([obj.params.numBeads,numel(n)], k(:,1),k(:,2));
            
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
            if exist('domainHandler','var')
           
                % The bead positions
               
                obj.positions.beads.prev.x(1) = 0;% the first bead is placed at the origin
                obj.positions.beads.prev.y(1) = 0;% the first bead is placed at the origin
                obj.positions.beads.prev.z(1) = 0;% the first bead is placed at the origin
                
                for bIdx = 2:obj.params.numBeads
                    inDomain= true;
                     while inDomain 
                         x = obj.params.b*randc(obj.params.dimension,1);
%                          n = sqrt(sum(x.^2));
                         for dIdx = 1:obj.params.dimension
                          obj.positions.beads.prev.(obj.dimNames{dIdx})(bIdx) = obj.positions.beads.prev.(obj.dimNames{dIdx})(bIdx-1)+x(dIdx);
                         end

                     a.x = obj.positions.beads.prev.x(bIdx);
                     a.y = obj.positions.beads.prev.y(bIdx);
                     a.z = obj.positions.beads.prev.z(bIdx);
                     inDomain = ~domainHandler.InDomain(a);
                     end                
               end
            
            else
                
            for dIdx = 1:obj.params.dimension
                % The bead positions
                obj.positions.beads.prev.(obj.dimNames{dIdx})(1) = 0;
                for bIdx = 2:obj.params.numBeads
                 obj.positions.beads.prev.(obj.dimNames{dIdx})(bIdx) = ...
                    obj.positions.beads.prev.(obj.dimNames{dIdx})(bIdx-1)+...
                    randn(1)*((obj.params.b)^(1/obj.params.dimension));
                end                
            end
            end
%             obj.positions.beads.prev = obj.positions.beads.cur;
        end                        
                
        function Next(obj,varargin)
            % Next simulation step 
%                 tic
                % Update simulation time and step 
%                 obj.time = obj.time + obj.params.dt;
                obj.step = obj.step + 1;
                obj.GetNewBeadsPosition; 
%                 obj.GetSprings; % set previous spring position as the current spring position
%                 obj.GetForces;  % calculate he forces in the current position
%                 endRoundTime = toc;
%                 obj.simulationTime.total        = obj.simulationTime.total+endRoundTime;
%                 obj.simulationTime.meanStepTime = (obj.simulationTime.meanStepTime*(obj.step-1)+endRoundTime)/obj.step;% the mean round time 
        end
        
        function GetForces(obj)% unused
             % Calculate forces   
             % call forceManager
             obj.GetSpringForce();             
             obj.GetLennardJonesForce();
             obj.GetBendingElasticityForce();
             obj.GetNoise();
        end        
        
        function Run(obj,varargin)
           
            while obj.runSimulation
               
                obj.Next                
                obj.runSimulation = obj.runSimulation & obj.step<obj.params.numSteps;
                
                continue
            end
        end
                
        function GetNewBeadsPosition(obj)
            % Solve the rouse system, given the previous location of the chain system
            % parameters
            % solve the system of equations
%             notify(obj,'getForce');
            
            prevPos      = [obj.positions.beads.prev.x,obj.positions.beads.prev.y,obj.positions.beads.prev.z]; 
%             ljForce      = [obj.forces.lennardJones.x, obj.forces.lennardJones.y, obj.forces.lennardJones.z];
%             bendingForce = [obj.forces.bending.x, obj.forces.bending.y,obj.forces.bending.z];
            
%             noiseForce   = [obj.forces.noise.x, obj.forces.noise.y,obj.forces.noise.z];
            p = obj.params;
            
            % Apply forces on the beads to get the new bead position
            newPos       = obj.handles.classes.forceManager.Apply(prevPos,obj.connectionMap.map,p.dt,...
                           p.diffusionConst,p.springConst,p.LJPotentialWidth,p.LJPotentialDepth,p.bendingConst,...
                           p.minBeadDist,p.fixedBeadNum);                             

             obj.positions.beads.cur.x = newPos(:,1);
             obj.positions.beads.cur.y = newPos(:,2);
             obj.positions.beads.cur.z = newPos(:,3);            
        end
        
        function Step(obj,force)% experimental 
            prevPos      = [obj.positions.beads.prev.x,obj.positions.beads.prev.y,obj.positions.beads.prev.z]; 
            ljForce      = force.lenardJones; %[obj.forces.lennardJones.x, obj.forces.lennardJones.y, obj.forces.lennardJones.z];
            bendingForce = force.bending;     % [obj.forces.bending.x, obj.forces.bending.y,obj.forces.bending.z];            
            noiseForce   = force.noise;       % [obj.forces.noise.x, obj.forces.noise.y,obj.forces.noise.z];
            

            newPos  = -obj.params.springConst.*obj.forces.springs*obj.params.dt*prevPos+...
                       ljForce+...
                       obj.params.bendingConst*obj.params.dt*bendingForce+...
                       noiseForce+...
                       prevPos;

             obj.positions.beads.cur.x = newPos(:,1);
             obj.positions.beads.cur.y = newPos(:,2);
             obj.positions.beads.cur.z = newPos(:,3);    
        end
        
        function SetPrevBeadPosition(obj)
             % Save previous position as the currentPosition                     
            obj.positions.beads.prev= obj.positions.beads.cur;
              
        end                                                                  
        
        function Reset(obj,newParams)% unfinished
            if nargin<2
                
            end
            
        end
        
    end
        
end

