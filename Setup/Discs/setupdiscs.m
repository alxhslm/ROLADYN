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
    D{i} = setup_each_disc(D{i},N,x{i});
    D{i}.iLocal = [(D{i}.iNode-1)*4 + (1:4) iDofCount + (1:D{i}.NDof)];
    iDofCount = iDofCount + D{i}.NDof;
end

function D = setup_each_disc(D,N,x0)

if ~isfield(D,'material')
    error('Missing field material');
end

D.z   = N(D.iNode);

D.material = setupmaterial(D.material);

if isfield(D,'m') 
    if isnan(D.material.rho)
        if isfield(D,'R') && isfield(D,'t')
            %only know geometric properties
            D = dim2mass(D);
        else
            error('Disc not fully defined');
        end
    elseif ~(isfield(D,'Id') && isfield(D,'R'))
        %only know inertia properties
        D = mass2dim(D);   
    end
elseif isnan(D.material.rho) && ~isfield(D,'R')
    error('Disc not fully defined');
end

required_fields = {'R','t'};
for i = 1:length(required_fields)
    if ~isfield(D,required_fields{i})
        error('Missing field %s from Disc',required_fields{i});
    end
end

if ~isfield(D,'bGyro')
    D.bGyro = 1;
end

D.NSegments = length(D.t);

if ~isfield(D,'NrSegment')
    D.NrSegment = 1;
end
if length(D.NrSegment) == 1
    D.NrSegment = D.NrSegment*ones(1,D.NSegments);
end

if ~isfield(D,'Nt')
    D.Nt = 4;
end

required_fields = {'material'};
for i = 1:length(required_fields)
    if ~isfield(D,required_fields{i})
        error('Missing field %s from Disc',required_fields{i});
    end
end

D.r = D.R(1);
D.iSegment = [];
for i = 1:D.NSegments
    [D.mSegment(i),D.IdSegment(i),D.IpSegment(i)] = disc_properties(D.material.rho,D.R(i),D.R(i+1),D.t(i));
    rSeg = linspace(D.R(i),D.R(i+1),D.NrSegment(i)+1);
    D.r = [D.r rSeg(2:end)];
    D.iSegment = [D.iSegment i*ones(1,D.NrSegment(i))];
end
D.dr = diff(D.r);
D.Nr = length(D.r);

D.mDisc = sum(D.mSegment);
D.IdDisc = sum(D.IdSegment);
D.IpDisc = sum(D.IpSegment);

D.theta = linspace(0,2*pi,D.Nt+1);
D.beta = diff(D.theta);
D.theta = D.theta(1:end-1);

params_required = {'iNode'};
for i = 1:length(params_required)
    if ~isfield(D,params_required{i})
        error('Cannot find parameter "%s" in the P.Rotor.Disc structure',params_required{i});
    end
end

%finite elements
NDofe = 3;
NEle = 4;
NHub = 4;
NRoot = 4;
D.NDof = D.Nt * D.Nr * NDofe + NHub;
NDofTot = D.NDof + NRoot;

for i = 1:D.Nt
    for j = 1:D.Nr
      D.SNode{i,j} = zeros(NDofe,NDofTot);
      D.SNode{i,j}(:,NHub+NRoot+((i-1)*D.Nr   + j - 1)*NDofe+(1:NDofe)) = eye(NDofe);
        
      %rotation matrix from Hub -> Disc Node coords
      D.RHub{i,j} = [0        0    D.r(j)*sin(D.theta(i))    -D.r(j)*cos(D.theta(i)); %w
                     0        0    D.r(j)*cos(D.theta(i))     D.r(j)*sin(D.theta(i))  %dw_dt
                     0        0     sin(D.theta(i))           -cos(D.theta(i))];      %dw_dr
    end
end
for i = 1:D.Nt
    if i < D.Nt
        ip1 = i+1;
    else
        ip1 = 1;
    end
    for j = 1:(D.Nr-1)
        A = D.beta(i)/2*((D.r(j)+D.dr(j))^2 - D.r(j)^2);
        iSeg = D.iSegment(j);
        
        D.me(i,j) = D.material.rho*A*D.t(iSeg);
        
        D.Se{i,j} = [D.SNode{i,j};
                     D.SNode{i,j+1};
                     D.SNode{ip1,j+1};
                     D.SNode{ip1,j}];
        
        D.Re{i,j} = eye(NDofe*NEle);
        
        [D.Ke{i,j},D.Me{i,j},D.Ge{i,j}] = annular(D.material,D.t(iSeg),D.theta(i),D.r(j),D.beta(i),D.dr(j));
        if ~D.bGyro
            D.Ge{i,j} = 0*D.Ge{i,j};
        end
        
        if ~isempty(x0)
            D.Fe0{i,j} = D.Ke{i,j}*D.Re{i,j}*D.Se{i,j}*x0;
        else
            D.Fe0{i,j} = zeros(NDofe*NEle,1);
        end
    end   
end    

%% Edge stiffness and damping
required_fields = {'KEdge_tt','KEdge_rr','KEdge_zz'};
for i = 1:length(required_fields)
    if ~isfield(D,required_fields{i})
        D.(required_fields{i}) = Inf;
    end
end
required_fields = {'CEdge_tt','CEdge_rr','CEdge_zz','KEdge_rt','CEdge_rt'};
for i = 1:length(required_fields)
    if ~isfield(D,required_fields{i})
        D.(required_fields{i}) = 0;
    end
end

D.KEdge = blkdiag(D.KEdge_zz, [D.KEdge_tt D.KEdge_rt;
                             D.KEdge_rt D.KEdge_rr]);
% D.KEdge = max(min(D.KEdge,1E20),-1E20);
ii = isinf(D.KEdge); D.KEdge(ii) = 0;

D.CEdge = blkdiag(D.CEdge_zz, [D.CEdge_tt D.CEdge_rt;
                      D.CEdge_rt D.CEdge_rr]);
D.CEdge = max(min(D.CEdge,1E20),-1E20);
ii = isinf(D.CEdge); D.CEdge(ii) = 0;

%% Hub inertia
D.SHub = zeros(NHub,NDofTot); D.SHub(:,NRoot+(1:NHub)) = eye(NHub);

required_fields = {'mHub','IdHub'};
for i = 1:length(required_fields)
    if ~isfield(D,required_fields{i})
        D.(required_fields{i}) = 0;
    end
end

if ~isfield(D,'IpHub')
    D.IpHub = 2*D.IdHub;
end

D.m  = D.mHub  + D.mDisc;
D.Id = D.IdHub + D.IdDisc;
D.Ip = D.IpHub + D.IpDisc;

D.MHub = diag([D.m D.m D.IdHub D.IdHub]);
D.GHub = D.bGyro*blkdiag(zeros(2),D.IpHub*antidiag([1 -1]));

%% Root compliance
D.SRoot = zeros(NRoot,NDofTot); D.SRoot(:,1:NRoot) = eye(NRoot);
required_fields = {'KRoot_rr','KRoot_tt'};
for i = 1:length(required_fields)
    if ~isfield(D,required_fields{i})
        D.(required_fields{i}) = Inf;
    end
end
default_fields = {'KRoot_rt','CRoot_rr','CRoot_rt','CRoot_tt'};
for i = 1:length(default_fields)
    if ~isfield(D,default_fields{i})
        D.(default_fields{i}) = 0;
    end
end
    
D.KRoot = [diag(D.KRoot_rr*[1 1])       D.KRoot_rt*eye(2);
              D.KRoot_rt*eye(2)  diag(D.KRoot_tt*[1 1])];
     
ii = isinf(D.KRoot); D.KRoot(ii) = 0;
% D.KRoot = max(min(D.KRoot,1E20),-1E20);

D.CRoot = [diag(D.CRoot_rr*[1 1])       D.CRoot_rt*eye(2);
              D.CRoot_rt*eye(2)  diag(D.CRoot_tt*[1 1])];
         
ii = isinf(D.CRoot); D.CRoot(ii) = 0;
% D.CRoot = max(min(D.CRoot,1E20),-1E20);

if ~isempty(x0)
    D.F0Root = D.KRoot*(D.SHub-D.SRoot)*x0;
else
    D.F0Root = zeros(4,1);
end

%% Totals

%compute the various matricies in the principal axes of the disk
D.M = diag([D.m D.m D.Id D.Id]);

%gyroscopic terms
D.G = D.bGyro*blkdiag(zeros(2),D.Ip*antidiag([1 -1]));