function R = setuprotors(R)
% Rotor mass and gyro matrices

if ~iscell(R)
    R = {R};
end

for i = 1:length(R)
    if ~isfield(R{i}, 'Name')
        R{i}.Name = sprintf('Rotor %d',i);
    end
    R{i} = setup_each_rotor(R{i});
end

function R = setup_each_rotor(R)

if ~isfield(R,'z')
    R.z = mean(R.Nodes);
end

%throw error if we don't have certain parameters
params_required = {'Nodes','Speed'};
for i = 1:length(params_required)
    if ~isfield(R,params_required{i})
        error('Cannot find parameter "%s" in the P.Rotor structure',params_required{i});
    end
end

%put a default rigid shaft, linking the discs/bearings
%if none specified

if ~isfield(R,'Shaft')
    R.Shaft = {};
end

if ~isfield(R,'Disc')
    R.Disc = {};
end

m  = 0;
Ip = 0;
mz = 0;

R.NDof = length(R.Nodes)*4;

%setup all of the shafts on the specified rotor
R.Shaft = setupshafts(R.Shaft, R.Nodes);
for i = 1:length(R.Shaft)
    m  = m  + R.Shaft{i}.Inertia.m;
    Ip = Ip + R.Shaft{i}.Inertia.Ip;
    mz = mz + R.Shaft{i}.Inertia.m * R.Shaft{i}.Inertia.z;
end

%and now the discs
R.Disc = setupdiscs(R.Disc,R.Nodes);
for i = 1:length(R.Disc)
    m  = m  + R.Disc{i}.Inertia.m;
    Ip = Ip + R.Disc{i}.Inertia.Ip;
    mz = mz + R.Disc{i}.Inertia.m * R.Disc{i}.z;
    R.NDof = R.NDof + R.Disc{i}.NDof;
end

R.Inertia.m = m;
R.Inertia.Ip = Ip;
R.Inertia.z = mz/m;

Id = 0;
for i = 1:length(R.Shaft)
    Id = Id + R.Shaft{i}.Inertia.Id + R.Shaft{i}.Inertia.m * (R.Shaft{i}.Inertia.z - R.Inertia.z)^2;
end
R.Inertia.Id = Id;

IMap = eye(R.NDof);
for i = 1:length(R.Nodes)
    R.SNode{i} = IMap((i-1)*4 + (1:4),:);
end