function [Params,F0,K,C,M] = setupRadial(Params)

if ~isfield(Params, 'c')
    Params.c = 0;
end

if ~isfield(Params, 'n')
    Params.n = 1;
end

if ~isfield(Params, 'Crad')
    Params.Crad = 0;
end

if ~isfield(Params, 'Krad')
    error('The stiffness field "Krad" is missing from the parameter structure')
end

if ~isfield(Params, 'Mrad')
    Params.Mrad = 0;
end

if ~isfield(Params, 'Crad')
    Params.Crad = 0;
end

if ~isfield(Params, 'Krot')
    error('The stiffness field "Krot" is missing from the parameter structure')
end

if ~isfield(Params, 'Mrot')
    Params.Mrad = 0;
end

F0 = zeros(8,1);
K = kron([1 -1; -1 1], diag([Params.Krad Params.Krot Params.Krad Params.Krot]));
C = kron([1 -1; -1 1], diag([Params.Crad Params.Crot Params.Crad Params.Crot]));
M = kron([1 -1; -1 1], diag([Params.Mrad Params.Mrot Params.Mrad Params.Mrot]));

Params.bActive = abs([Params.Krad Params.Krot Params.Krad Params.Krot]) > 0;
Params.bRigid  = isinf([Params.Krad Params.Krot Params.Krad Params.Krot]);

% Clip Infs
Params.Krad = max(-1E20,min(Params.Krad,1E20));
Params.Krot = max(-1E20,min(Params.Krot,1E20));
