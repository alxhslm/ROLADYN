function [Params,F0,K,C,M] = setupEmpty(Params)
F0 = zeros(8,1);
K = zeros(8);
C = zeros(8);
M = zeros(8);
Params.bActive = false(4,1);
Params.bRigid = false(4,1);