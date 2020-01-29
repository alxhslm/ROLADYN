function D = setupdiscs(D,N)
if ~iscell(D)
    D = {D};
end

iDofCount = length(N)*4;
for i = 1:length(D)
    if ~isfield(D{i}, 'Name')
        D{i}.Name = sprintf('Disc %d',i);
    end
    D{i} = setup_each_disc(D{i},N);
    D{i}.iLocal = [(D{i}.iNode-1)*4 + (1:4) iDofCount + (1:D{i}.NDof)];
    iDofCount = iDofCount + D{i}.NDof;
end

function D = setup_each_disc(D,N)

if ~isfield(D,'iNode')
    error('Missing field "iNode" from disc "%s"',D.Name);
end

D.z = N(D.iNode);

if ~isfield(D,'Type')
    error('Missing field "Type" from disc "%s"',D.Name);
end

if ~isfield(D,'Ring')
    error('Missing field "Ring from disc "%s"',D.Name);
end

%% Options
if ~isfield(D,'Options')
    D.Options = struct();
end
if ~isfield(D.Options,'bGyro')
    D.Options.bGyro = 1;
end

%ndof
D.Hub.NDof = 4;
D.Root.NDof = 4;

%% Ring
%actual body of disc
switch D.Type
    case 'Rigid'
        D = setupRigidDisc(D);
    case 'Flexible'
        D = setupFlexibleDisc(D);
end

%% Hub inertia
%central mass to which the rings connect
D.Hub.S = zeros(D.Hub.NDof,D.NDofTot); D.Hub.S(:,D.Root.NDof+(1:D.Hub.NDof)) = eye(D.Hub.NDof);

required_fields = {'m','Id','Ip'};
for i = 1:length(required_fields)
    if ~isfield(D.Hub,required_fields{i})
        D.Hub.(required_fields{i}) = 0;
    end
end

D.Hub.M = diag([D.Hub.m D.Hub.m D.Hub.Id D.Hub.Id]);
D.Hub.G = D.Options.bGyro*blkdiag(zeros(2),D.Hub.Ip*antidiag([1 -1]));

%% Root compliance
%flexible joint between hub and shaft
D.Root.S = zeros(D.Root.NDof,D.NDofTot); D.Root.S(:,1:D.Root.NDof) = eye(D.Root.NDof);
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

%% Totals
D.Inertia.m  = D.Hub.m  + sum(D.Ring.Inertia.m);
D.Inertia.Id = D.Hub.Id + sum(D.Ring.Inertia.Id);
D.Inertia.Ip = D.Hub.Ip + sum(D.Ring.Inertia.Ip);

%compute the various matrices in the principal axes of the disk
D.M = diag([D.Inertia.m D.Inertia.m D.Inertia.Id D.Inertia.Id]);

%gyroscopic terms
D.G = D.Options.bGyro*blkdiag(zeros(2),D.Inertia.Ip*antidiag([1 -1]));

function D = setupFlexibleDisc(D)

if ~isfield(D.Ring.Geometry,'t')
    error('Need disc thickness "t" for disc "%s" of type "Flexible"',D.Name);
end
if ~isfield(D.Ring.Geometry,'R')
    error('Need disc radius "R" for disc "%s" of type "Flexible"',D.Name);
end

if length(D.Ring.Geometry.R) == length(D.Ring.Geometry.t)
    D.Ring.Geometry.R = [0 D.Ring.Geometry.R];
end

if length(D.Ring.Geometry.R) ~= (length(D.Ring.Geometry.t) + 1)
    error('Mismatching lengths for "R" and "t" for disc "%s" of type "Flexible"',D.Name);
end

NSegments = length(D.Ring.Geometry.t);

%number of radial points per ring
if ~isfield(D.Ring,'Nr')
    D.Ring.Nr = 1;
end
if length(D.Ring.Nr) == 1
    D.Ring.Nr = D.Ring.Nr*ones(1,NSegments);
end

% Compute mass of each ring
D.Ring.m = zeros(NSegments,1);
D.Ring.Id = zeros(NSegments,1);
D.Ring.Ip = zeros(NSegments,1);
for i = 1:NSegments
    [D.Ring.m(i),D.Ring.Id(i),D.Ring.Ip(i)] = disc_properties(D.Material.rho,D.Ring.R(i),D.Ring.R(i+1),D.Ring.t(i));
end

%totals
D.Ring.Inertia.m  = sum(D.Ring.m);
D.Ring.Inertia.Id = sum(D.Ring.Id);
D.Ring.Inertia.Ip = sum(D.Ring.Ip);
        
%need to store radial inertia of each ring, as this is not accounted for by FE model
D.Ring.M = diag([D.Ring.Inertia.m  D.Ring.Inertia.m  0 0]);
D.Ring.G = zeros(4);

%% Finite element model
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
D.NDof = D.Mesh.Nt * D.Mesh.Nr * NDofe + D.Hub.NDof;
D.NDofTot = D.NDof + D.Root.NDof;

theta = linspace(0,2*pi,D.Mesh.Nt+1);
D.Mesh.dtheta = diff(theta);
D.Mesh.theta  = theta(1:end-1);

for i = 1:D.Mesh.Nt
    for j = 1:D.Mesh.Nr
      D.Mesh.SNode{i,j} = zeros(NDofe,D.NDofTot);
      D.Mesh.SNode{i,j}(:,D.Hub.NDof+D.Root.NDof+((i-1)*D.Mesh.Nr   + j - 1)*NDofe+(1:NDofe)) = eye(NDofe);

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

%% Edge stiffness
% connection to inner hub

%stiffness and damping
if ~isfield(D,'Edge')
    D.Edge = struct();
end

required_fields = {'Ktt','Krr','Kzz'};
for i = 1:length(required_fields)
    if ~isfield(D.Edge,required_fields{i})
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
ii = isinf(D.Edge.K); D.Edge.K(ii) = 0;

D.Edge.C = blkdiag(D.Edge.Czz, [D.Edge.Ctt D.Edge.Crt;
                                D.Edge.Crt D.Edge.Crr]);
ii = isinf(D.Edge.C); D.C(ii) = 0;

function D = setupRigidDisc(D)

if isfield(D.Ring,'Inertia')
    %have inertial properties -- just use them
    required_fields = {'m','Id','Ip'};
    for i = 1:length(required_fields)
        if ~isfield(D.Ring.Inertia,required_fields{i})
            error('Need disc inertia "%s" for disc "%s" of type "Rigid"',required_fields{i},D.Name);
        end
    end
    
    if ~isfield(D.Ring,'Geometry')
        D.Ring.Geometry = struct();
    end

    optional_fields = {'R','t'};
    default_val = {0.1,0.02};
    for i = 1:length(optional_fields)
        if ~isfield(D.Ring.Geometry,optional_fields{i})
            warning('Default disc diemsion "%s" for disc "%s" of type "Rigid"',optional_fields{i},D.Name);
            D.Ring.(optional_fields{i}) = default_val(i);
        end
    end
elseif  isfield(D.Ring,'Geometry')
    if ~isfield(D,'Material')
        error('Missing field "Material" from disc "%s" of type "%s"',D.Name,D.Type);
    end
    D.Material = setupmaterial(D.Material);
    D = setupRingGeometry(D);
else
    error('Not enough information given for disc "%s" of type "Rigid"',D.Name);
end

%need to store radial inertia of each ring, as this is not accounted for by FE model
D.Ring.M = diag([sum(D.Ring.Inertia.m) sum(D.Ring.Inertia.m) sum(D.Ring.Inertia.Id) sum(D.Ring.Inertia.Id)]);
D.Ring.G = D.Options.bGyro*blkdiag(zeros(2),sum(D.Ring.Inertia.Ip)*antidiag([1 -1]));

D.NDof = D.Hub.NDof;
D.NDofTot = D.NDof + D.Root.NDof;

function D = setupRingGeometry(D)
if ~isfield(D.Ring.Geometry,'t')
    error('Need disc thickness "t" for disc "%s" of type "Flexible"',D.Name);
end
if ~isfield(D.Ring.Geometry,'R')
    error('Need disc radius "R" for disc "%s" of type "Flexible"',D.Name);
end

if length(D.Ring.Geometry.R) == length(D.Ring.Geometry.t)
    D.Ring.Geometry.R = [0 D.Ring.Geometry.R];
end

if length(D.Ring.Geometry.R) ~= (length(D.Ring.Geometry.t) + 1)
    error('Mismatching lengths for "R" and "t" for disc "%s" of type "Flexible"',D.Name);
end

NSegments = length(D.Ring.Geometry.t);

% Compute mass of each ring
D.Ring.Inertia.m = zeros(NSegments,1);
D.Ring.Inertia.Id = zeros(NSegments,1);
D.Ring.Inertia.Ip = zeros(NSegments,1);
for i = 1:NSegments
    [D.Ring.Inertia.m(i),D.Ring.Inertia.Id(i),D.Ring.Inertia.Ip(i)] = disc_properties(D.Material.rho,D.Ring.Geometry.R(i),D.Ring.Geometry.R(i+1),D.Ring.Geometry.t(i));
end