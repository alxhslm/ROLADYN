function S = setupshafts(S,N)
if ~iscell(S)
    S = {S};
end
for i = 1:length(S)
    if ~isfield(S{i}, 'Name')
        S{i}.Name = sprintf('Shaft %d',i);
    end
    S{i} = setup_each_shaft(S{i},N);
    S{i}.iLocal = ((S{i}.iNodes(:)-1)*4 + (1:4))';
    S{i}.iLocal = S{i}.iLocal(:);
end

function S = setup_each_shaft(S,N)
if ~isfield(S, 'Name')
    S.Name = sprintf('Shaft %d',ind);
end

if ~isfield(S,'iNodes')
    error('Missing field "iNodes" from shaft "%s"',S.Name);
end

if ~isfield(S,'Material')
    error('Missing field "Material" from shaft "%s"',S.Name);
end
S.Material = setupmaterial(S.Material);

%% Section properties
if ~isfield(S,'Section')
    error('Missing field "Section" from shaft "%s"',S.Name);
end
if ~isfield(S.Section,'ri')
    S.Section.ri = 0;
end
S.Section.A = pi*(S.Section.ro^2 - S.Section.ri^2);
S.Section.I = pi/4*(S.Section.ro^4 - S.Section.ri^4);

%% Options
if ~isfield(S,'Options')
    S.Options = struct();
end
if ~isfield(S.Options,'bGyro')
    S.Options.bGyro = 1;
end
if ~isfield(S.Options,'Element')
    S.Options.Element = 'euler_bernoulli';
end

%% Finite element model
RShaft = eye(8);
RShaft = RShaft([1 4 5 8 2 3 6 7],:); %rotor to bearing coords
RShaft([6 8],:) = -1*RShaft([6 8],:); %flip signs of angles in xz plane

% Create 1D mesh
S.Mesh.z   = N(S.iNodes);
S.Mesh.Nz  = length(S.Mesh.z);
IMap = eye(4*length(S.iNodes));
for i = 1:length(S.iNodes)
    S.Mesh.SNode{i} = IMap((i-1)*4 + (1:4),:);
end

for k = 1:(S.Mesh.Nz-1)
    %element mass, stiffness etc
    S.Element{k}.L = S.Mesh.z(k+1) - S.Mesh.z(k);
    S.Element{k}.z = 0.5*(S.Mesh.z(k+1) + S.Mesh.z(k));
    S.Element{k}.Inertia.m = S.Material.rho*S.Section.A*S.Element{k}.L;
    S.Element{k}.Inertia.Id = S.Element{k}.Inertia.m/12*S.Element{k}.L^2 + S.Material.rho*S.Section.I*S.Element{k}.L;
    S.Element{k}.Inertia.Ip = 2*S.Material.rho*S.Section.I*S.Element{k}.L;
    
    [S.Element{k}.K,S.Element{k}.M,S.Element{k}.G] = feval(['shaft_' S.Options.Element],S.Material,S.Section.ri,S.Section.ro,S.Element{k}.L);
    S.Element{k}.G = S.Options.bGyro*S.Element{k}.G;
    
    S.Element{k}.R = RShaft;
    S.Element{k}.S = [S.Mesh.SNode{k}; S.Mesh.SNode{k+1}];
end

%% Totals
S.Section.L = 0;
S.Inertia.m = 0;
S.Inertia.Ip = 0;
mz = 0;
for k = 1:length(S.Element)
    S.Section.L = S.Section.L + S.Element{k}.L;
    S.Inertia.m = S.Inertia.m + S.Element{k}.Inertia.m;
    S.Inertia.Ip = S.Inertia.Ip + S.Element{k}.Inertia.Ip;
    mz = mz + S.Element{k}.Inertia.m * S.Element{k}.z;
end
S.Inertia.z = mz / (S.Inertia.m + eps);

S.Inertia.Id = 0;
for k = 1:length(S.Element)
    S.Inertia.Id = S.Inertia.Id + S.Element{k}.Inertia.Id + S.Element{k}.Inertia.m * (S.Element{k}.z - S.Inertia.z)^2;
end

S.M = diag([S.Inertia.m S.Inertia.m S.Inertia.Id S.Inertia.Id]);
S.G = S.Options.bGyro*blkdiag(zeros(2),antidiag(S.Inertia.Id*[1 -1]));