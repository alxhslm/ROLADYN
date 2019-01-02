function results = rotor_backbone(P,hbm,A,w0,iDofX,iDofU,x0)
%FRF_NONLIN Summary of this function goes here
%   Detailed explanation goes here

problem = rotor2problem(P);
P = problem.P;
hbm = setuphbm(hbm,problem);

problem.iDof = iDofX;
problem.iInput = iDofU;

if nargin < 7
    x0 = rotor_init(hbm,problem,w0,A(1));
end
% problem.Xscale = max(max(abs(x0)),1E-6);
results = hbm_bb(hbm,problem,A(1),A(end),w0,x0);
results.O = results.w;
results = rmfield(results,'w');
results.x0 = x0;