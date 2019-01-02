function U = rotor_hbm_excite(hbm,problem,w0)
P = problem.P;

Fe = excitation_hbm(P);
U = packharm3d(Fe,hbm.harm);

function F = packharm3d(Fh,harm)
kHarm = harm.kHarm;
kHarm(:,1) = kHarm(:,1) * harm.rFreqRatio(1) * harm.rFreqBase(1);
kHarm(:,2) = kHarm(:,2) * harm.rFreqRatio(2) * harm.rFreqBase(2);

F = zeros(size(kHarm,1),size(Fh{1},1));
i1 = find(kHarm(:,1) == 1 & kHarm(:,2) == 0);
if ~isempty(i1)
    F(i1,:) = Fh{1}.';
end
i2 = find(kHarm(:,1) == 0 & kHarm(:,2) == 1);
if ~isempty(i2)
    F(i2,:) = Fh{2}.';
end