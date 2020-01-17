function S = setupstators(S,x0)
% Rotor mass and gyro matrices

if ~iscell(S)
    S = {S};
end

for i = 1:length(S)
    if nargin > 1
        x{i} = S{i}.S*x0;
    else
        x{i} = [];
    end
    S{i} = setup_each_stator(S{i},i,x{i});
end

function S = setup_each_stator(S,ind,x0)

if ~isfield(S, 'Name')
    S.Name = sprintf('Rotor %d',ind);
end

if ~isfield(S,'NDof')
    S.NDof = size(S.M,1);
end

if ~isempty(x0)
    S.F0 = S.K*x0;
else
    S.F0 = zeros(S.NDof,1);
end

fields = {'M','K','C'};
for i = 1:length(fields)
    if ~isfield(S,fields{i})
        S.(fields{i}) = zeros(S.NDof);
    end
end

S.Ks = S.K;
ii = isinf(S.K); S.K(ii) = 0;
ii = isinf(S.C); S.C(ii) = 0;