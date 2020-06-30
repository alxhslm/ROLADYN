function [D,K,C,M] = setupSFD(D)
if ~isfield(D,'Nt')
    D.Nt = 11; 
end
if ~isfield(D,'Nz')
    D.Nz = 13; 
end
D.NDof = 0;

if ~isfield(D,'bPinned')
    D.bPinned = 0;
end

if size(D.KbSquirrel,1) == 2 
    D.KbSquirrel = [diag([D.KbSquirrel(1,1) Inf]) diag([D.KSbquirrel(1,2) 0]);
                    diag([D.KSbquirrel(2,1) 0])   diag([D.KbSquirrel(2,2) Inf])]; 
end

D.fun = str2func(['SFD_', D.Model]);

D.bRigid = isiinf(diag(D.KbSquirrel));

D.KSq = max(-1E20,min(D.KbSquirrel,1E20));
K = D.K0 + D.KbSquirrel;
C = D.C0;
M = D.M0;