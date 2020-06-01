function [Piezo,K,C,M] = setupPiezo(Piezo)
Piezo.Model = 'goldfarb';

Piezo.fun = str2func(['piezo_', Piezo.Model]);

Piezo.Mech.kO = Piezo.Mech.k + Piezo.Elec.T^2/Piezo.Elec.C;
Piezo.Mech.kS = Piezo.Mech.k + Piezo.Elec.T^2/(Piezo.Elec.C + Piezo.Elec.Cm);

if ~isfield(Piezo,'bAmplifiersOn')
    Piezo.bAmplifiersOn = true;
end

if Piezo.bAmplifiersOn
    Piezo.NDofTot = 2;
    kpz = Piezo.Mech.kS;
else
    Piezo.NDofTot = 0;
    kpz = Piezo.Mech.kO;
end

%assemble outputs
K = kron([1 -1; -1 1], diag([kpz,Inf,Inf,Inf]));
C = kron([1 -1; -1 1], diag([Piezo.Mech.c, 0,0,0]));
M = kron([1 -1; -1 1], zeros(4));