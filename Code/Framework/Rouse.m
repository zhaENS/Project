classdef Rouse<handle
    
    % TODO: sort the properties by subjects, move unneccessary to private
    %      properties 
    % TODO: add chain starting position options: linear, dense, loose, random,
    %       etc...
    % TODO: when beads are connected, make their stating position close.
    % make the starting position dense by restricting each bead to lay at a
    % distance from the center of mass (sequentially in creating the
    % initial beads)
    % TODO: add Reset option (so the class doesn't have to be initialized again)
    % TODO: add an option to recreate the noise values from previously used
    %       seed or simulation parameters
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
        distributions  
        connectionMap

        runSimulation
        eq        
        mobilityMatrices
        
        relaxationTime %[sec] relaxation time of the first mode
        
        beadsDist
        forces
        dimNames = {'x','y','z'};
    end
    
    properties (Access=private)
        noise
        noiseIndex = 0;
        numNoisePoints
    end
    
    events
    end
    
    methods
        function obj = Rouse(rouseParams)
            % class constructor
            obj.runSimulation = true;
            obj.params        = rouseParams;
          
%             obj.SetDefaultParams
            % adjust input parameters to match the structure expected in
            % the class
            obj.SetInputParams
            
            obj.CalculateRelaxationTime
            
            obj.InitializeRouseStruct
            
%             obj.CreateControls             
            
        end
        
        function SetDefaultParams(obj)%Unused
            % Default Params
            
            obj.params.bendingElasticityForce = true;
            obj.params.lennardJonesForce   = true;
            obj.params.springForce       = true;

            obj.params.dimension         = 3;       % an integer between 1 and 3
            obj.params.numBeads          = 2;       % number of beads in the chain            
            obj.params.dt                = 10^(-2); % time step [sec]
            obj.params.k                 = 1e-1;    % Boltzman constant            
            obj.params.b                 = 1;            
            obj.params.zeta              = 1;             
            obj.params.temperature       = 273;      % thermodynamic temperature (Kelvin)            
            obj.params.minBeadDist       = 0.5;      % min distance between beads [unused in this version]            
            obj.params.fixedBeadNum      = [];       % beads index fixed            
            obj.params.noiseDistribution = 'Gaussian';            
            obj.params.diffusionConst    = 1;            
            obj.params.springConst       = -2*obj.params.dimension*obj.params.diffusionConst*obj.params.dt/(obj.params.b^2);            
            obj.params.noiseStd          = sqrt(2*obj.params.diffusionConst*obj.params.dt);%2*obj.params.k*obj.params.temperature*obj.params.dt;            
            obj.params.noiseMean         = 0;            
            obj.params.LJPotentialDepth  = 1e-5;             
            obj.params.LJPotentialWidth  = 10;            
            obj.params.bendingConst      = 0.1;                        
            obj.params.connectMonomers   = [];% bounded monomers should come in pairs, this is an integer list of [n by 2]
            obj.params.beta              = 2*ones(obj.params.numBeads); % the anomolouse exponent. for beta~=2 all beads are connected. 
            
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
            
            obj.positions.beads.cm.x            = 0;
            obj.positions.beads.cm.y            = 0;
            obj.positions.beads.cm.z            = 0;
            
            obj.positions.springs.x             = sparse(1:obj.params.numBeads,1:obj.params.numBeads,0); 
            obj.positions.springs.y             = sparse(1:obj.params.numBeads,1:obj.params.numBeads,0);  
            obj.positions.springs.z             = sparse(1:obj.params.numBeads,1:obj.params.numBeads,0); 
            obj.positions.springs.lengths       = sparse(1:obj.params.numBeads,1:obj.params.numBeads,0);
            
            % The angle between springs is defined as a 3D array
            % but is represented as N 2D arrays of sparse matrices 
            obj.positions.springs.angleBetweenSprings = sparse(repmat(1:obj.params.numBeads,1,obj.params.numBeads),1:obj.params.numBeads^2,pi);
            obj.beadsDist                       = sparse(1:obj.params.numBeads,1:obj.params.numBeads,0);
            obj.distributions.beadDistance.dist = []; % to be moved to distributionHandler
            obj.distributions.beadDistance.mean = []; % remove
            obj.distributions.beadDistance.std  = []; % remove
            obj.distributions.beadDistance.rms  = []; % remove
            obj.distributions.endToEnd.vec      = []; % remove
            obj.distributions.endToEnd.std      = []; % remove
            obj.distributions.endToEnd.mean     = []; % remove
            obj.distributions.endToEnd.dist     = []; % remove
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
            obj.SetBeadsMobilityMatrices;  
            obj.SetBeadConnectionMap
            obj.SetInitialChainPosition; 
            obj.SetBeadsEquation 
        end
        
        function SetBeadsEquation(obj)% unused
            % the beads positions 
%             p       - previous position
%             Fm      - forward mobility matrix
%             Bm      - backward mobility matrix 
%             Fb      - forward binary matrix
%             Bb      - backward binary matrix
%             n       - noise term 
%             b       - mean springs length
% %             params  - parameters 
            if obj.params.minBeadDist~=0

                obj.eq.beads = @(p,Fb,Bb,mobilityMatrices,D,b,dt,n,l,dim)...
                    (-dim*D*dt/(b^2))*abs(mobilityMatrices.beads.forward*p)+...
                    (-dim*D*dt/(b^2))*abs(mobilityMatrices.beads.backward*(p))-...
                    (-dim*2*D*dt/(b^2))*l*ones(size(p))+...
                    n+...
                    p;
            else
                % for the case the min distance is 0, the position can be
                % found solving the matrix form 
                obj.eq.beads = @(prevPos,mobilityMatrix,D,b,dt,n,dim)(-dim*D*dt/(b.^2))*(mobilityMatrix.beads.normal*prevPos)+n+prevPos;
            end
                
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
            obj.SetRouseMatrix                         
            
            [v,lambda]           = obj.GetRouseEig;
            % if the lambdas are changed, then the connection map is
            % changed. 
            % adjust connect monomer list accordingly 
            obj.AdjustConnectedMonomersList
            cm = obj.params.connectMonomers;
            assert(all(cm(:)<=obj.params.numBeads),'The monomers index to connect cannot exceed the number of monomers in the chain')
           
            % recalculate the connection map 
            connectionMatrix = v*diag(lambda)*v';
            connectionMatrix = connectionMatrix.*(abs(connectionMatrix)>1e-12);
            for bIdx = 1:size(cm,1)
                connectionMatrix(cm(bIdx,1),cm(bIdx,2))  = -1;
                connectionMatrix(cm(bIdx,2),cm(bIdx,1))  = -1;
            end
            
            obj.connectionMap.connectionMatrix = connectionMatrix; % used for multiplication of beads location 
            cMap = double((connectionMatrix-diag(diag(connectionMatrix)))~=0);
            cMap(~cMap) = NaN;
           % add the connection between distant beads
           
           for bIdx = 1:size(cm,1)
               cMap(cm(bIdx,1),cm(bIdx,2))= 1;
               cMap(cm(bIdx,2),cm(bIdx,1))= 1;
           end
           
           obj.connectionMap.map                = cMap;
           [obj.connectionMap.indices.in.map]   = ~isnan(cMap); % linear indices of positions of monomer pairs connected. 
           [obj.connectionMap.indices.in.list(:,1), obj.connectionMap.indices.in.list(:,2)] = find(triu(~isnan(cMap)));
           [obj.connectionMap.indices.out.map]  =  isnan(cMap);% linear indices of positions of monomer pairs connected. 
           [obj.connectionMap.indices.out.list(:,1), obj.connectionMap.indices.out.list(:,2)]  = find(obj.connectionMap.indices.out.map); 
           cMapNz = cMap;
           cMapNz(isnan(cMapNz))=0;
           obj.connectionMap.numConnections = sum(cMapNz,2);% number of connection each bead has 

           % Set angle between springs as a triplet of beads. The middle
           % bead represents the middle node in the triplet defining the
           % angle. The convention is row number, column number, height
           % number. To represent an angle between bead i,j,k. start from
           % row i to column j and go to the depth layer k
           
           sMap = zeros(obj.params.numBeads,obj.params.numBeads,obj.params.numBeads);
           % first register the linear connections 
           for bIdx = 2:obj.params.numBeads-1
               sMap(bIdx-1,bIdx,bIdx+1) = 1;
               obj.connectionMap.beadTriplets(end+1,:) = [bIdx-1,bIdx,bIdx+1];
           end
           
           % next, register the non-trivial connections
%            cMap1 = triu(cMap);
%            cMap1(isnan(cMap1))=0;
%            cMap1(isnan(cMap1))=0;
           for bIdx = 1:obj.params.numBeads
               % there are either 2 or more non-zero elements in each row.
               allPairs = [];
                f= find(cMap(:,bIdx)==1);
                for fIdx = 1:numel(f);
                    g = find(cMap(f(fIdx),:)==1);
                    g= setdiff(g,bIdx);
                    for gIdx = 1:numel(g)
                        allPairs(end+1,:) = [min([bIdx g(gIdx)]),f(fIdx),max([bIdx,g(gIdx)])];
                    end
                end

%                 if numel(f)>2
%                     allPairs = [f,bIdx*ones(numel(f),1),circshift(f,1)];
%                 elseif numel(f)==2
%                     allPairs = [f(1),bIdx,f(2)];                    
%                 else 
%                     allPairs = [];
%                 end
                % sort the indices such that the smaller index is on the
                % left
%                 for pIdx = 1:size(allPairs,1)
%                     if allPairs(pIdx,1)>allPairs(pIdx,3)
%                         allPairs(pIdx,:) = fliplr(allPairs(pIdx,:));
%                     end                        
%                 end
%                 allPairs = unique(allPairs,'rows');
                for pIdx = 1:size(allPairs,1)
                    % make a symmetric matrix in 3D
                    sMap(allPairs(pIdx,1),allPairs(pIdx,2),allPairs(pIdx,3))= 1;
                    sMap(allPairs(pIdx,3),allPairs(pIdx,2), allPairs(pIdx,1))=1;% check                     
                    obj.connectionMap.beadTriplets(end+1,:) = allPairs(pIdx,:);
                end
              obj.connectionMap.angles = sMap;  
           end
            obj.connectionMap.beadTriplets= unique(obj.connectionMap.beadTriplets,'rows');
           % define the linear indices in the array of bead triplets to
           % save running time later in the calculation of angle between
           % beads. angles is the angle at the middle bead (bead2) so we
           % have bead1ToBead2 and bead3ToBead2 to define the edges of the
           % angle
           bt     = obj.connectionMap.beadTriplets;
           for btIdx = 1:numel(bt(:,1))
            obj.connectionMap.indices.in.linear.bead1ToBead2(btIdx,1)  = sub2ind([obj.params.numBeads,obj.params.numBeads],min([bt(btIdx,1),bt(btIdx,2)]),max([bt(btIdx,1),bt(btIdx,2)]));
            obj.connectionMap.indices.in.linear.bead3ToBead2(btIdx,1)  = sub2ind([obj.params.numBeads,obj.params.numBeads],min([bt(btIdx,2),bt(btIdx,3)]),max([bt(btIdx,2),bt(btIdx,3)]));
           end
           % bead triplets represents the indices in the N 2D sparse matrix
           % of bead connections, row for bead1, bead2 is represented N
           % times in N*2D matrix, bead3 is the column Nth block meaning
           % size(matrix,2)/numBeads.                       
           obj.connectionMap.indices.in.linear.beadTriplets  = (obj.connectionMap.beadTriplets(:,3)-1)*obj.params.numBeads^2+obj.connectionMap.indices.in.linear.bead1ToBead2;
           k = sub2ind([obj.params.numBeads, obj.params.numBeads],obj.connectionMap.indices.in.list(:,1),obj.connectionMap.indices.in.list(:,2));
           obj.connectionMap.indices.in.linear.beadPairs  = k;% pair of beads connected (linear indices of the matrix size [numBead X numBeads])
           [c,~,ic] = unique(obj.connectionMap.beadTriplets(:,2),'R2012a');                
           obj.connectionMap.indices.uniqueBeadTriplets.uniqueList = c;
           obj.connectionMap.indices.uniqueBeadTriplets.uniqueListIndices = ic;
           n       = 1:numel(obj.connectionMap.beadTriplets(:,2));
           k       = [c(ic),n'];  
           obj.connectionMap.indices.uniqueBeadTriplets.uniqueIndicesMap = k;
           n = [1:obj.params.numBeads,1:abs(numel(obj.connectionMap.beadTriplets(:,2))-obj.params.numBeads)];
           k       = obj.connectionMap.indices.uniqueBeadTriplets.uniqueIndicesMap;
           obj.connectionMap.indices.uniqueBeadTripletsLinear = sub2ind([obj.params.numBeads,numel(n)], k(:,1),k(:,2));
            
        end
        
        function SetRouseMatrix(obj)
            % Define the Rouse matrix 
            obj.mobilityMatrices.rouse          = 2*eye(obj.params.numBeads)...
                                                  -diag(ones(1,obj.params.numBeads-1),1)...
                                                  -diag(ones(1,obj.params.numBeads-1),-1);
            obj.mobilityMatrices.rouse(1,1)     = 1;
            obj.mobilityMatrices.rouse(end,end) = 1;
        end
        
        function AdjustConnectedMonomersList(obj)% unfinished
            % This function operates after the Rouse Matrix is defined. It
            % takes into account the connected monomers indices defined by
            % the user and the resulted connection matrix after the beta
            % values of the chains have been altered (Reminder: if beta=2
            % we get the classical Rouse matrix, other wise it is the beta
            % model) 
            m= obj.mobilityMatrices.rouse;
            % remove the diagonal 
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
            % Initialize the springs vectors and lengths 
            obj.GetSprings;
%             % calculate forces at initial position 
            obj.GetForces;
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
        
        function GetForces(obj)
             % Calculate forces               
             obj.GetSpringForce();             
             obj.GetLennardJonesForce();
             obj.GetBendingElasticityForce();
             obj.GetNoise();
        end
        
        function SetBeadsMobilityMatrices(obj)%unused
            % create the Rouse mobility matrix for the beads
            % this matrix is used in the calculation of the new beads positions
            numBeads              = obj.params.numBeads;
            
            forwardMobilityMatrix = diag(ones(1,numBeads),0) +...
                diag(-1*ones(1,numBeads-1),1);                
            forwardMobilityMatrix(end,:) = 0;            
            obj.mobilityMatrices.beads.forward = forwardMobilityMatrix;
            
            % set backward mobility matrix 
            backwardMobilityMatrix = diag(ones(1,numBeads),0) +...
                diag(-1*ones(1,numBeads-1),-1);                
            backwardMobilityMatrix(1,:) = 0;            
            obj.mobilityMatrices.beads.backward = backwardMobilityMatrix;
            
            % set normal (minDist =0 ) mobility matrix 
            obj.mobilityMatrices.beads.normal = 2*diag(ones(1,numBeads),0) -diag(ones(1,numBeads-1),-1)-...
                diag(ones(1,numBeads-1),+1);
            obj.mobilityMatrices.beads.normal(1,1)     = 1;
            obj.mobilityMatrices.beads.normal(end,end) = 1;
            
            if ~isempty(obj.params.fixedBeadNum)
                f = (obj.params.fixedBeadNum); %#ok<*ST2NM>
                for fIdx = 1:numel(f)
                    obj.mobilityMatrices.beads.normal(f(fIdx),:) = 0;                    
                end
            end
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
            prevPos      = [obj.positions.beads.prev.x,obj.positions.beads.prev.y,obj.positions.beads.prev.z]; 
            ljForce      = [obj.forces.lennardJones.x, obj.forces.lennardJones.y, obj.forces.lennardJones.z];
            bendingForce = [obj.forces.bending.x, obj.forces.bending.y,obj.forces.bending.z];
            
            noiseForce   = [obj.forces.noise.x, obj.forces.noise.y,obj.forces.noise.z];
            
            newPos  = -(obj.params.diffusionConst/obj.params.b^2)*obj.params.dimension*obj.params.dt*obj.forces.springs*prevPos+...
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
        
        function GetSpringForce(obj)
           % get the force between beeds caused by spring
           % attrcation/repultion  
            obj.GetSprings;
%             obj.forces.springs = zeros(size(l));
            if obj.params.springForce
                L                  = 1-obj.params.minBeadDist./obj.positions.springs.lengths;
                L(obj.connectionMap.indices.out.map)= 0; % zero out forces where there is no connection
                sumForces          = sum(L,2);
                force              = diag(sumForces)-L; % set the maindiagonal
                obj.forces.springs = force;
                obj.forces.springs(obj.params.fixedBeadNum,:) = 0;
                obj.forces.springs(:,obj.params.fixedBeadNum) = 0;
%                 obj.forces.springs = (obj.forces.springs~=0).*obj.connectionMap.connectionMatrix;
            end
        end
        
        function GetLennardJonesForce(obj)
            % Get the LJ forces between monomers
                        
            if obj.params.lennardJonesForce
                beadsDistMat = (obj.beadsDist);                               
                sig          = obj.params.LJPotentialWidth;
                epsilon      = obj.params.LJPotentialDepth;
                bdmInv       = (beadsDistMat).^(-1);            % one over the bead distance matrix 
                d            = obj.MatIntPower(bdmInv, 6); % matrix integer power 
                t            = (sig^6).*d;                
                forceValue   = 24*(epsilon*bdmInv).*(-2*t.*t +t); % derivative of LJ function 
                
                forceValue(isnan(forceValue))= 0;
                for dIdx = 1:3
                    % replicate the position vector
                    A    = obj.positions.beads.cur.(obj.dimNames{dIdx});
                    siz  = [1,obj.params.numBeads];
                    B1   = A(:,ones(siz(2),1));
                    siz  = [obj.params.numBeads,1];
                    A    = A';
                    B2   = A(ones(siz(1),1),:);
                    % Subtract positions to get the direction vectors
                    fd  = B1-B2;

                     obj.forces.lennardJones.(obj.dimNames{dIdx}) = (sum(fd.*forceValue,1)')*obj.params.dt;
                     obj.forces.lennardJones.(obj.dimNames{dIdx})(obj.params.fixedBeadNum) = 0; 
                end
            end
            
        end
        
        function GetBendingElasticityForce(obj)
            % Calculate the force and direction 

            if (obj.params.bendingElasticityForce)
                obj.CalculateAngleBetweenSprings
                % get the coordinates of the first vector comprising the
                % angle
                t1x     = obj.positions.springs.x(obj.connectionMap.indices.in.linear.bead1ToBead2);
                t1y     = obj.positions.springs.y(obj.connectionMap.indices.in.linear.bead1ToBead2);
                t1z     = obj.positions.springs.z(obj.connectionMap.indices.in.linear.bead1ToBead2);
                % get the coordinates of the second vector  comprising the
                % angle
                t2x     = obj.positions.springs.x(obj.connectionMap.indices.in.linear.bead3ToBead2);
                t2y     = obj.positions.springs.y(obj.connectionMap.indices.in.linear.bead3ToBead2);
                t2z     = obj.positions.springs.z(obj.connectionMap.indices.in.linear.bead3ToBead2);
               
                % calculate the the force (normalized) given by the angle
                FValue = (pi- obj.positions.springs.angleBetweenSprings(obj.connectionMap.indices.in.linear.beadTriplets))/pi;
                % direct the force of each bead toward the base of the
                % triangle formed by the two edges (vectors) comprising the
                % angle
                tempFx  = (t1x-(0.5*(t2x-t1x))).*FValue;
                tempFy  = (t1y-(0.5*(t2y-t1y))).*FValue;
                tempFz  = (t1z-(0.5*(t2z-t1z))).*FValue;  
%                 % ====== debug================
%                  hold on 
%                  line('XData',obj.positions.beads.cur.x(2),...
%                       'YData',obj.positions.beads.cur.y(2),...
%                       'ZData',obj.positions.beads.cur.z(2),...
%                       'Marker','o',...
%                       'MarkerFaceColor','r',...
%                       'Tag','point')
%                  % plot the vectors 
%                  quiver3(obj.positions.beads.cur.x(2),...
%                      obj.positions.beads.cur.y(2),...
%                      obj.positions.beads.cur.z(2),t1x,t1y,t1z)
%                  quiver3(obj.positions.beads.cur.x(2),...
%                      obj.positions.beads.cur.y(2),...
%                      obj.positions.beads.cur.z(2),t2x,t2y,t2z)
%                   quiver3(obj.positions.beads.cur.x(2),...
%                      obj.positions.beads.cur.y(2),...
%                      obj.positions.beads.cur.z(2),tempFx,tempFy,tempFz)
%                  hold off
%                   
%                  %=== end debug ====================
                 
                % for each bead, calculate the force as the sum of forces
                % exerted by the angles
%                 n       = [1:obj.params.numBeads,1:abs(numel(obj.connectionMap.beadTriplets(:,2))-obj.params.numBeads)];
%                 a       = sparse(obj.params.numBeads,1:numel(n),0);     
                a       = sparse(size(obj.positions.springs.angleBetweenSprings,1),size(obj.positions.springs.angleBetweenSprings,2));
                ind     = obj.connectionMap.indices.uniqueBeadTripletsLinear;
                a       = full(a);
                a(ind)  = tempFx;

                Force.x = sum(a,2);
                Force.x(obj.params.fixedBeadNum) = 0;
                a(ind)  = tempFy;
                Force.y = sum(a,2);
                Force.y(obj.params.fixedBeadNum) = 0;
                a(ind)  = tempFz;
                Force.z = sum(a,2);
                Force.z(obj.params.fixedBeadNum) = 0;
                obj.forces.bending = Force;               
            end
           
        end        
        
        function GetNoise(obj)
            % Get noise according to the specified noise districution params
            % the output noise term is a colum vector wuth params.numBeads entries  
%             noiseTerm = zeros(obj.params.numBeads,3);
            if strcmpi(obj.params.noiseDistribution,'Gaussian')
                
                if obj.step==1 || mod(obj.step,1e5)==0;
                    obj.noise = obj.params.noiseMean+...
                    randn(1e5*obj.params.numBeads,obj.params.dimension)*obj.params.noiseStd;
                    obj.noiseIndex = 1;
                else 
                    obj.noiseIndex = obj.noiseIndex+1;
                end
%                 noiseTerm(:,1:obj.params.dimension) = obj.params.noiseMean+...
%                     randn(obj.params.numBeads,obj.params.dimension)*obj.params.noiseStd;
                noiseTerm(:,1:obj.params.dimension) = obj.noise((obj.noiseIndex-1)*obj.params.numBeads +1: obj.noiseIndex*obj.params.numBeads,:);
                noiseTerm((obj.params.fixedBeadNum),:) = 0;        
               
                obj.forces.noise.x = noiseTerm(:,1);
                obj.forces.noise.y = noiseTerm(:,2);
                obj.forces.noise.z = noiseTerm(:,3);
            else
                error('Noise distribution is undefined')
            end
        end
              
        function CalculateAngleBetweenSprings(obj)
            % Calculate the angle between springs
            
            linInd1 = obj.connectionMap.indices.in.linear.bead1ToBead2;
            linInd2 = obj.connectionMap.indices.in.linear.bead3ToBead2;
            c1      = full([obj.positions.springs.x(linInd1),obj.positions.springs.y(linInd1),obj.positions.springs.z(linInd1)]);
            c2      = full([obj.positions.springs.x(linInd2),obj.positions.springs.y(linInd2),obj.positions.springs.z(linInd2)]);
%             c1 = c1(sum(c1,2)~=0,:);
%             c2 = c2(sum(c2,2)~=0,:);
            sAngle  = acos(sum(c1.*c2,2)./(sqrt(sum(c1.^2,2)).*sqrt(sum(c2.^2,2))));
            linInd3 = obj.connectionMap.indices.in.linear.beadTriplets;%(obj.connectionMap.beadTriplets(:,3)-1)*obj.params.numBeads^2+linInd1;
            obj.positions.springs.angleBetweenSprings(linInd3) = sAngle;

        end
                       
        function GetChainCenterOfMass(obj)
            obj.positions.beads.cm.x = mean(obj.positions.beads.cur.x(:,end));
            obj.positions.beads.cm.y = mean(obj.positions.beads.cur.y(:,end));
            obj.positions.beads.cm.z = mean(obj.positions.beads.cur.z(:,end));
        end
        
        function GetSprings(obj)
            % Calculate the springs vectors. the springs vectors are in the
            % direction of the bead_N-bead_N-1
            if obj.params.numBeads>1   

                for dIdx = 1:3
                 obj.positions.springs.(obj.dimNames{dIdx})(obj.connectionMap.indices.in.linear.beadPairs) = ...
                     obj.positions.beads.prev.(obj.dimNames{dIdx})(obj.connectionMap.indices.in.list(:,1))-...
                     obj.positions.beads.prev.(obj.dimNames{dIdx})(obj.connectionMap.indices.in.list(:,2));
                end
                
                % Get pairwise bead distances (between all beads)
                 obj.GetBeadsDist
                                                   
                 obj.positions.springs.lengths = obj.beadsDist.*obj.connectionMap.map;
            else
                obj.positions.springs.x       = 0;
                obj.positions.springs.y       = 0;
                obj.positions.springs.z       = 0;
                obj.positions.springs.lengths = 0;
            end                        
        end    
        
        function GetBeadsDist(obj)
           % Calculate the pairwise distance between beads
           % This function requires the mex file pdist2mex in the working path

            obj.beadsDist = pdist2mex([obj.positions.beads.prev.x,obj.positions.beads.prev.y, obj.positions.beads.prev.z]',...
                [obj.positions.beads.prev.x,obj.positions.beads.prev.y, obj.positions.beads.prev.z]','euc',[],[],[]);
            
        end
        
        function Reset(obj,newParams)% unfinished
            if nargin<2
                
            end
            
        end
        
    end
    
    methods (Static)

            function d = MatIntPower(a, n)
                % d = MatIntPower(a, n)
                % power with integer exponent
                % return the array d = a.^n

                % binary decomposition
                nb = dec2bin(n)-'0';
                ak = a;
                if nb(end)
                    d = a;
                else
                    d = 1;
                end
                for nbk = nb(end-1:-1:1)
                    ak = ak.*ak;
                    if nbk
                        d = d.*ak;
                    end
                end

            end % MatIntPower
       end
    
end

