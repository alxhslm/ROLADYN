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

if ~isfield(S,'NDof')
    S.NDof = size(S.M,1);
end

if ~isfield(S, 'Name')
    S.Name = sprintf('Rotor %d',ind);
end

if ~isempty(x0)
    S.F0 = S.K*x0;
else
    S.F0 = zeros(S.NDof,1);
end

S.Ks = S.K;

ii = isinf(S.K); S.K(ii) = 100;
ii = isinf(S.C); S.C(ii) = 100;