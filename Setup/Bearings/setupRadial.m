function [B,K,C,M] = setupRadial(B)

if ~isfield(B, 'c')
    B.c = 0;
end

if ~isfield(B, 'n')
    B.n = 1;
end

if ~isfield(B, 'C')
    B.C = 0;
end

if ~isfield(B, 'M')
    B.M = 0;
end

if ~isfield(B, 'K')
    error('The stiffness field "K" is missing from the parameter structure')
end

if ~isfield(B, 'KbParallel')
    B.KbParallel = zeros(4);
end

B.bActive = [true; false; true; false];
B.bActive = B.bActive | diag(B.KbParallel) > 0;

B.bRigid = isinf(diag(B.KbParallel));

B.KPar = kron([1 -1; -1 1], max(-1E20,min(B.KbParallel,1E20)));
K = kron([1 -1; -1 1], B.KbParallel);
C = kron([1 -1; -1 1],B.C*diag([1 0 1 0]));
M = kron([1 -1; -1 1],B.M*diag([1 0 1 0]));