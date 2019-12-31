function problem = rotor2problem(P)
problem.model = @rotor_hbm_model;
problem.excite = @rotor_hbm_excite;

if ~isfield(P.Model,'bCompressREB')
    P.Model.bCompressREB = 0;
end
P.Model.Reduced = get_internal_indices(P.Bearing,P.Model.bCompressREB);

NDofInt = P.Model.Reduced.NDofInt;
problem.sparsity = [ones(P.Model.NDof)    ones(P.Model.NDof,NDofInt);
                    ones(NDofInt,P.Model.NDof) P.Model.Reduced.Sparsity];

problem.NDof    = P.Model.NDof;
problem.NOutput = P.Model.NDof;
problem.NInput  = P.Model.Excite.NExcite;

problem.NDof = problem.NDof + NDofInt;

problem.Ku = P.Model.Excite.K;
problem.Cu = P.Model.Excite.C;
problem.Mu = P.Model.Excite.M;

problem.M  =  P.Model.Rotor.M + P.Model.Stator.M;
problem.C  =  P.Model.Rotor.C + P.Model.Stator.C;
problem.K  =  P.Model.Rotor.K + P.Model.Stator.K;
problem.F0 = -P.Model.Fg;

problem.Ku = [problem.Ku; zeros(NDofInt,size(problem.Ku,2))];
problem.Cu = [problem.Cu; zeros(NDofInt,size(problem.Cu,2))];
problem.Mu = [problem.Mu; zeros(NDofInt,size(problem.Mu,2))];
problem.K  = blkdiag(problem.K,zeros(NDofInt));
problem.C  = blkdiag(problem.C,zeros(NDofInt));
problem.M  = blkdiag(problem.M,zeros(NDofInt));
problem.F0 = [problem.F0; zeros(NDofInt,1)];

if isfield(P.Model,'RDofPlot')
    problem.RDofPlot = P.Model.RDofPlot;
elseif isfield(P.Model,'iDofPlot')
    problem.iDofPlot = P.Model.iDofPlot;
end

if ~isfield(P.Model,'bUseGroups')
    P.Model.bUseGroups = 0;
end

if P.Model.bUseGroups
    problem.iGroup = [1*ones(P.Model.NDof,1);
                      2*ones(NDofInt,1)];
else
    problem.iGroup = [1*ones(P.Model.NDof,1);
                      1*ones(NDofInt,1)];
end

if ~isfield(P.Model,'bAnalyticalDerivs')
    P.Model.bAnalyticalDerivs = 1;
end

problem.P = P;

function Red = get_internal_indices(B,bCompressREB)
iDofIn = 0;
iInt = [];
iSum = {};
for i = 1:length(B)
    if strcmp(B{i}.Model,'REB')
        REB = B{i}.Params;
        if bCompressREB
            iInt = [iInt; iDofIn+((1:REB.Model.NDof)-1)*REB.Setup.Z+1];
            for k = 1:REB.Model.NDof
                iSum = [iSum; iDofIn+ (k-1)*REB.Setup.Z + (1:REB.Setup.Z)];
            end
            iDofIn = iDofIn + B{i}.NDofInt;
            Sparsity{i} = ones(B{i}.Params.Model.NDof);
        else
            iInt = [iInt; iDofIn+(1:B{i}.NDofInt)];
            iSum = [iSum; num2cell(iDofIn+(1:B{i}.NDofInt))];
            iDofIn  = iDofIn  + B{i}.NDofInt;
            Sparsity{i} = repmat(eye(REB.Setup.Z),B{i}.Params.Model.NDof);
        end
    else
        %othwerise need all the internal states
        iInt = [iInt; iDofIn+(1:B{i}.NDofInt)];
        iSum = [iSum; num2cell(iDofIn+(1:B{i}.NDofInt)')];
        iDofIn  = iDofIn  + B{i}.NDofInt;
        Sparsity{i} = ones(B{i}.NDofInt);
    end
    
end
Red.NDofInt = length(iInt);
Red.iInt = iInt;
Red.iSum = iSum;

Sparsity = Sparsity';
Red.Sparsity = blkdiag(Sparsity{:});