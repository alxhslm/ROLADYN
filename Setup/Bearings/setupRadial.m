function [B,F,K,C,xInt] = setupRadial(B,qi,qo,Oi,Oo)

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

States.qi = qi;
States.qo = qo;
States.qidot = 0*qi;
States.qodot = 0*qo;
States.bSolve = 1;
[Forces,~,Stiffness] = radial_model(B, States);
B.F0 = Forces.F;
B.K0 = Stiffness.K;
B.C0 = Stiffness.C;
B.qi0 = qi;
B.qo0 = qo;

if ~isfield(B, 'KbParallel')
    B.KbParallel = zeros(4);
end

B.bActive = [true; false; true; false];
B.bActive = B.bActive | diag(B.KbParallel) > 0;

B.KPar = kron([1 -1; -1 1], max(-1E20,min(B.KbParallel,1E20)));
K = B.K0 + kron([1 -1; -1 1], B.KbParallel);
C = B.C0;
F = B.F0 + B.KPar*[qi; qo];

xInt = [];