function [B,K,C] = setupRadial(B)

if ~isfield(B, 'c')
    B.c = 0;
end

if ~isfield(B, 'n')
    B.n = 1;
end

if ~isfield(B, 'C')
    B.C = 0;
end

if ~isfield(B, 'K')
    error('The stiffness field "K" is missing from the parameter structure')
end

if ~isfield(B, 'KbParallel')
    B.KbParallel = zeros(4);
end

B.bActive = [true; false; true; false];
B.bActive = B.bActive | diag(B.KbParallel) > 0;

B.KPar = kron([1 -1; -1 1], max(-1E20,min(B.KbParallel,1E20)));
K = B.K0 + kron([1 -1; -1 1], B.KbParallel);
C = B.C0;