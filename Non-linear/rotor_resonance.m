function [q,f,u,w] = resonance(P,hbm,A,w0,iDofX,iDofU,q0)
%FRF_NONLIN Summary of this function goes here
%   Detailed explanation goes here

problem.model = @rotor_hbm;
problem.NOutput = 4*length(P.Bearing)*2;
problem.NInput = P.Mesh.NExcite;
problem.NDof   = P.Model.NDofTot;

problem.P = P;
problem.Ku = [P.Model.Excite.Kub*P.Mesh.Excite.Sub; zeros(P.Model.NDofInt,P.Mesh.NExcite)];
problem.Cu = [P.Model.Excite.Cub*P.Mesh.Excite.Sub; zeros(P.Model.NDofInt,P.Mesh.NExcite)];
problem.Mu = [P.Model.Excite.Mub*P.Mesh.Excite.Sub; zeros(P.Model.NDofInt,P.Mesh.NExcite)];
problem.M  = blkdiag(P.Model.Rotor.M,zeros(P.Model.NDofInt));
problem.C  = blkdiag(P.Model.Rotor.C,zeros(P.Model.NDofInt));
problem.K  = blkdiag(P.Model.Rotor.K,zeros(P.Model.NDofInt));
problem.F0 = [-P.Model.Fg + P.Model.Rotor.F0 - P.Model.Rotor.K*P.Model.x0; zeros(P.Model.NDofInt,1)];

if strcmp(P.Bearing{1}.Model{1},'REB')
    REB = P.Bearing{1}.Params{1};
    problem.sparsity = [ones(P.Model.NDof) ones(P.Model.NDof,REB.Z*REB.NDof);
                        ones(REB.Z*REB.NDof,P.Model.NDof) repmat(eye(REB.Z),REB.NDof)];
else
    problem.sparsity = ones(P.Model.NDof);
end

problem.iDof = iDofX;
problem.iInput = iDofU;

if nargin < 7
    x0 = rotor_init(hbm,problem,w0,A);
else
    x0 = q0 * P.Model.A;
end
problem.Xscale = max(max(abs(x0)),1E-6);

[X,w,u,f] = hbm_resonance(hbm,problem,w0,A,x0);
q = X(:,1:P.Model.NDof)*P.Model.A';
