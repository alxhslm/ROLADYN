function [Params,F0,K,C,M] = setupSFD(Params)
if ~isfield(Params,'Nt')
    Params.Nt = 11; 
end
if ~isfield(Params,'Nz')
    Params.Nz = 13; 
end
Params.NDof = 0;

Params.fun = str2func(['SFD_', Params.Model]);

Params.bActive = true(4,1);
Params.bRigid = false(4,1);

F0 = zeros(8,1);
K = zeros(8);
C = zeros(8);
M = zeros(8);