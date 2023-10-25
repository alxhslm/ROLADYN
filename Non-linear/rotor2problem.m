function problem = rotor2problem(P)
problem.model = @rotor_hbm_model;
problem.excite = @rotor_hbm_excite;

if ~isfield(P.Model,'bCompressREB')
    P.Model.bCompressREB = 0;
end

problem.sparsity = ones(P.Model.NDof);

problem.NDof    = P.Model.NDof;
problem.NOutput = P.Model.NDof;
problem.NInput  = P.Model.Excite.NExcite;

problem.Ku = P.Model.Excite.K;
problem.Cu = P.Model.Excite.C;
problem.Mu = P.Model.Excite.M;

problem.M  =  P.Model.Rotor.M + P.Model.Stator.M + P.Model.Bearing.Lin.M;
problem.C  =  P.Model.Rotor.C + P.Model.Stator.C + P.Model.Bearing.Lin.C;
problem.K  =  P.Model.Rotor.K + P.Model.Stator.K + P.Model.Bearing.Lin.K;
problem.G  =  P.Model.Rotor.G;
problem.F0 = -P.Model.Fg + P.Model.Bearing.Lin.F;

if isfield(P.Model,'RDofPlot')
    problem.RDofPlot = P.Model.RDofPlot;
elseif isfield(P.Model,'iDofPlot')
    problem.iDofPlot = P.Model.iDofPlot;
end

if ~isfield(P.Model,'bAnalyticalDerivs')
    P.Model.bAnalyticalDerivs = 1;
end

problem.P = P;
