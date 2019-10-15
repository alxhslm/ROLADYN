function D = setupdiscs(D,N,x0)
if ~iscell(D)
    D = {D};
end

iDofCount = length(N)*4;
for i = 1:length(D)
    if  nargin > 2 && ~isempty(x0)
        x{i} = D{i}.S*x0;
    else
        x{i} = [];
    end
    if ~isfield(D{i}, 'Name')
        D{i}.Name = sprintf('Disc %d',i);
    end
    D{i} = setup_each_disc(D{i},N,x{i});
    D{i}.iLocal = [(D{i}.iNode-1)*4 + (1:4) iDofCount + (1:D{i}.NDof)];
    iDofCount = iDofCount + D{i}.NDof;
end

function D = setup_each_disc(D,N,x0)

if ~isfield(D,'iNode')
    error('Missing field "iNode" from disc "%s"',D.Name);
end

D.z = N(D.iNode);

if ~isfield(D,'Material')
    error('Missing field "Material" from disc "%s"',D.Name);
end
D.Material = setupmaterial(D.Material);

if ~isfield(D,'Ring')
    D.Ring = struct();
end
if ~isfield(D.Ring,'R')
    D.Ring.R = NaN(1,2);
end
if ~isfield(D.Ring,'t')
    D.Ring.t = NaN;
end

if ~isfield(D,'Inertia')
    D.Inertia = struct();
end

if any(isnan(D.Ring.R)) || any(isnan(D.Ring.t))
    if length(D.Ring.t) > 1
        error('Unable to compute disc geometry for disc "%s"',D.Name);
    elseif isfield(D,'Inertia')
        D = mass2dim(D);
    else
        error('Need either field "Inertia" or "Ring" for disc "%s"',D.Name);
    end
end

%% Options
if ~isfield(D,'Options')
    D.Options = struct();
end
if ~isfield(D.Options,'bGyro')
    D.Options.bGyro = 1;
end

%% Compute mass of each ring
if ~isfield(D,'Ring')
    D.Ring = struct();
end
if ~isfield(D.Ring,'t')
    D.Ring.t = [];
end
NSegments = length(D.Ring.t);

%number of radial points per ring
if ~isfield(D.Ring,'Nr')
    D.Ring.Nr = 1;
end
if length(D.Ring.Nr) == 1
    D.Ring.Nr = D.Ring.Nr*ones(1,NSegments);
end

D.Ring.m = zeros(NSegments,1);
D.Ring.Id = zeros(NSegments,1);
D.Ring.Ip = zeros(NSegments,1);
for i = 1:NSegments
    [D.Ring.m(i),D.Ring.Id(i),D.Ring.Ip(i)] = disc_properties(D.Material.rho,D.Ring.R(i),D.Ring.R(i+1),D.Ring.t(i));
end

%need to store radial inertia of each ring, as this is not accounted for by FE model
D.Ring.M = diag([sum(D.Ring.m) sum(D.Ring.m) 0 0]);

%% Finite element model
NHub = 4;
NRoot = 4;
if NSegments > 0
    %create mesh
    D.Mesh.r = D.Ring.R(1);
    D.Mesh.t = [];

    for i = 1:NSegments
        rSeg = linspace(D.Ring.R(i),D.Ring.R(i+1),D.Ring.Nr(i)+1);
        tSeg = ones(1,D.Ring.Nr(i))*D.Ring.t(i);
        D.Mesh.r = [D.Mesh.r rSeg(2:end)];
        D.Mesh.t = [D.Mesh.t tSeg];
    end
    D.Mesh.dr = diff(D.Mesh.r);
    D.Mesh.Nr = length(D.Mesh.r);

    if ~isfield(D.Mesh,'Nt')
        D.Mesh.Nt = 4;
    end

    % transformation matrices
    NDofe = 3;
    NEle = 4;
    D.NDof = D.Mesh.Nt * D.Mesh.Nr * NDofe + NHub;
    NDofTot = D.NDof + NRoot;

    theta = linspace(0,2*pi,D.Mesh.Nt+1);
    D.Mesh.dtheta = diff(theta);
    D.Mesh.theta  = theta(1:end-1);

    for i = 1:D.Mesh.Nt
        for j = 1:D.Mesh.Nr
          D.Mesh.SNode{i,j} = zeros(NDofe,NDofTot);
          D.Mesh.SNode{i,j}(:,NHub+NRoot+((i-1)*D.Mesh.Nr   + j - 1)*NDofe+(1:NDofe)) = eye(NDofe);

          %rotation matrix from Hub -> Disc Node coords
          D.Mesh.RHub{i,j} = [0        0    D.Mesh.r(j)*sin(D.Mesh.theta(i))    -D.Mesh.r(j)*cos(D.Mesh.theta(i)); %w
                              0        0    D.Mesh.r(j)*cos(D.Mesh.theta(i))     D.Mesh.r(j)*sin(D.Mesh.theta(i))  %dw_dt
                              0        0     sin(D.Mesh.theta(i))                -cos(D.Mesh.theta(i))];      %dw_dr
        end
    end
    
    % elements mass etc
    for i = 1:D.Mesh.Nt
        iNext = mod(i,D.Mesh.Nt)+1;
        for j = 1:(D.Mesh.Nr-1)
            A = D.Mesh.dtheta(i)/2*(D.Mesh.r(j+1)^2 - D.Mesh.r(j)^2);

            D.Element{i,j}.m = D.Material.rho*A*D.Mesh.t(j);

            D.Element{i,j}.S = [D.Mesh.SNode{i,j};
                                D.Mesh.SNode{i,j+1};
                                D.Mesh.SNode{iNext,j+1};
                                D.Mesh.SNode{iNext,j}];

            D.Element{i,j}.R = eye(NDofe*NEle);

            [D.Element{i,j}.K,D.Element{i,j}.M,D.Element{i,j}.G] = disc_annular(D.Material,D.Mesh.t(j),D.Mesh.r(j),D.Mesh.dtheta(i),D.Mesh.dr(j));

            D.Element{i,j}.G = D.Options.bGyro*D.Element{i,j}.G;

            if ~isempty(x0)
                D.Element{i,j}.F0 = D.Element{i,j}.K*D.Element{i,j}.R*D.Element{i,j}.S*x0;
            else
                D.Element{i,j}.F0 = zeros(NDofe*NEle,1);
            end
        end   
    end    
else
    D.Mesh.Nt = [];
    D.Mesh.Nr = [];
    D.NDof = NHub;
    NDofTot = D.NDof + NRoot;
end

%% Hub
%stiffness and damping
if ~isfield(D,'Hub')
    D.Hub = struct();
end

required_fields = {'Ktt','Krr','Kzz'};
for i = 1:length(required_fields)
    if ~isfield(D.Hub,required_fields{i})
        D.Edge.(required_fields{i}) = Inf;
    end
end
required_fields = {'Ctt','Crr','Czz','Krt','Crt'};
for i = 1:length(required_fields)
    if ~isfield(D.Edge,required_fields{i})
        D.Edge.(required_fields{i}) = 0;
    end
end

D.Edge.K = blkdiag(D.Edge.Kzz, [D.Edge.Ktt D.Edge.Krt;
                              D.Edge.Krt D.Edge.Krr]);
% D.Edge.K = max(min(D.Edge.K,1E20),-1E20);
ii = isinf(D.Edge.K); D.Edge.K(ii) = 0;

D.Edge.C = blkdiag(D.Edge.Czz, [D.Edge.Ctt D.Edge.Crt;
                                D.Edge.Crt D.Edge.Crr]);
D.Edge.C = max(min(D.Edge.C,1E20),-1E20);
ii = isinf(D.Edge.C); D.C(ii) = 0;

% inertia
D.Hub.S = zeros(NHub,NDofTot); D.Hub.S(:,NRoot+(1:NHub)) = eye(NHub);

required_fields = {'m','Id','Ip'};
for i = 1:length(required_fields)
    if ~isfield(D.Hub,required_fields{i})
        D.Hub.(required_fields{i}) = 0;
    end
end

D.Hub.M = diag([D.Hub.m D.Hub.m D.Hub.Id D.Hub.Id]);
D.Hub.G = D.Options.bGyro*blkdiag(zeros(2),D.Hub.Ip*antidiag([1 -1]));

%% Root compliance
D.Root.S = zeros(NRoot,NDofTot); D.Root.S(:,1:NRoot) = eye(NRoot);
required_fields = {'Krr','Ktt'};
for i = 1:length(required_fields)
    if ~isfield(D.Root,required_fields{i})
        D.Root.(required_fields{i}) = Inf;
    end
end
default_fields = {'Krt','Crr','Crt','Ctt'};
for i = 1:length(default_fields)
    if ~isfield(D.Root,default_fields{i})
        D.Root.(default_fields{i}) = 0;
    end
end
    
D.Root.K = [diag(D.Root.Krr*[1 1])   D.Root.Krt*eye(2);
              D.Root.Krt*eye(2)    diag(D.Root.Ktt*[1 1])];
     
ii = isinf(D.Root.K); D.Root.K(ii) = 0;
% D.Root.K = max(min(D.Root.K,1E20),-1E20);

D.Root.C = [diag(D.Root.Crr*[1 1])   D.Root.Crt*eye(2);
              D.Root.Crt*eye(2)    diag(D.Root.Ctt*[1 1])];
         
ii = isinf(D.Root.C); D.Root.C(ii) = 0;
% D.Root.C = max(min(D.Root.C,1E20),-1E20);

if ~isempty(x0)
    D.Root.F0 = D.Root.K*(D.Hub.S-D.Root.S)*x0;
else
    D.Root.F0 = zeros(4,1);
end

%% Totals
D.Inertia.m  = D.Hub.m  + sum(D.Ring.m);
D.Inertia.Id = D.Hub.Id + sum(D.Ring.Id);
D.Inertia.Ip = D.Hub.Ip + sum(D.Ring.Ip);

%compute the various matricies in the principal axes of the disk
D.M = diag([D.Inertia.m D.Inertia.m D.Inertia.Id D.Inertia.Id]);

%gyroscopic terms
D.G = D.Options.bGyro*blkdiag(zeros(2),D.Inertia.Ip*antidiag([1 -1]));