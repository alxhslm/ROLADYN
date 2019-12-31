function [Piezo,K,C,M] = setupPiezo(Piezo)
Piezo.Model = 'goldfarb';

Piezo.fun = str2func(['piezo_', Piezo.Model]);

Piezo.Mech.kO = Piezo.Mech.k + Piezo.T^2/Piezo.C;
Piezo.Mech.kS = Piezo.Mech.k + Piezo.T^2/(Piezo.C + Piezo.Cm);

%assemble outputs
K = kron([1 -1; -1 1], diag([Piezo.Mech.kS,Inf,Inf,Inf]));
C = kron([1 -1; -1 1], diag([Piezo.Mech.c, Inf,Inf,Inf]));
M = kron([1 -1; -1 1], diag([0, Inf,Inf,Inf]));

Piezo.NDofTot = 1;