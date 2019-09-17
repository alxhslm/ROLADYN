function R = setuprotors(R,x0)
% Rotor mass and gyro matrices

if ~iscell(R)
    R = {R};
end

for i = 1:length(R)
    if nargin > 1
        x{i} = R{i}.S*x0;
    else
        x{i} = [];
    end
    R{i} = setup_each_rotor(R{i},i,x{i});
end

function R = setup_each_rotor(R,ind,x0)

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
    R.Shaft = struct();
end

if ~isfield(R, 'Name')
    R.Name = sprintf('Rotor %d',ind);
end

m  = 0;
Id = 0;
Ip = 0;
mz = 0;

R.NDof = length(R.Nodes)*4;

%setup all of the shafts on the specified rotor
R.Shaft = setupshafts(R.Shaft, R.Nodes, x0);
for i = 1:length(R.Shaft)
    m  = m  + R.Shaft{i}.m;
    Id = Id + R.Shaft{i}.Id;
    Ip = Ip + R.Shaft{i}.Ip;
    mz = mz + R.Shaft{i}.m * mean(R.Shaft{i}.z);
end

%and now the discs
if isfield(R,'Disc')
    R.Disc = setupdiscs(R.Disc, R.Nodes, x0);
    for i = 1:length(R.Disc)
        m  = m  + R.Disc{i}.Inertia.m;
        Id = Id + R.Disc{i}.Inertia.Id;
        Ip = Ip + R.Disc{i}.Inertia.Ip;
        mz = mz + R.Disc{i}.Inertia.m * R.Disc{i}.z;
        R.NDof = R.NDof + R.Disc{i}.NDof;
    end
else
    R.Disc = {};
end

R.m = m;
R.Id = Id;
R.Ip = Ip;
R.z = mz/m;

IMap = eye(R.NDof);
for i = 1:length(R.Nodes)
    R.SNode{i} = IMap((i-1)*4 + (1:4),:);
end