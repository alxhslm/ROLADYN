function problem = rotor2problem(P)
problem.model = @rotor_hbm_model;
problem.excite = @rotor_hbm_excite;

P.Model.Reduced = get_internal_indices(P.Bearing,P.Model.bCompressREB);

NDofInt = P.Model.Reduced.NDofInt;
problem.sparsity = [ones(P.Model.NDof)    ones(P.Model.NDof,NDofInt);
                    ones(NDofInt,P.Model.NDof) P.Model.Reduced.Sparsity];

problem.NDof = P.Model.NDof;
problem.NOutput = P.Model.NDof;
problem.NInput = P.Mesh.NExcite;

if P.Model.bUseAlgebraic
    problem.NAlg   = NDofInt;
else
    problem.NAlg = 0;
    problem.NDof = problem.NDof + NDofInt;
end

problem.res.iDof = 1;
problem.res.iHarm = 2;

if P.Model.bCompressREB
    problem.jacobX = @rotor_hbm_jacobian;
%     problem.derivW = @rotor_hbm_deriv;
end

problem.P = P;

problem.Ku = P.Model.Excite.Kub*P.Mesh.Excite.Sub;
problem.Cu = P.Model.Excite.Cub*P.Mesh.Excite.Sub;
problem.Mu = P.Model.Excite.Mub*P.Mesh.Excite.Sub;

problem.M  = P.Model.Rotor.M;
problem.C  = P.Model.Rotor.C + P.Model.Bearing.C;
problem.K  = P.Model.Rotor.K;
problem.F0 = -P.Model.Fg + P.Model.Rotor.F0 - P.Model.Rotor.K*P.Model.x0;

if ~P.Model.bUseAlgebraic
    problem.Ku = [problem.Ku; zeros(NDofInt,size(problem.Ku,2))];
    problem.Cu = [problem.Cu; zeros(NDofInt,size(problem.Cu,2))];
    problem.Mu = [problem.Mu; zeros(NDofInt,size(problem.Mu,2))];
    problem.K  = blkdiag(problem.K,zeros(NDofInt));
    problem.C  = blkdiag(problem.C,zeros(NDofInt));
    problem.M  = blkdiag(problem.M,zeros(NDofInt));
    problem.F0 = [problem.F0; zeros(NDofInt,1)];
end

if isfield(P.Model,'iDofPlot')
    problem.iDofPlot = P.Model.iDofPlot;
else
    problem.iDofPlot = 1:P.Model.NDofTot;
end

if P.Model.bUseGroups
    problem.iGroup = [1*ones(P.Model.NDof,1);
                      2*ones(NDofInt,1)];
else
    problem.iGroup = [1*ones(P.Model.NDof,1);
                      1*ones(NDofInt,1)];
end

% problem.Ku = [P.Model.Excite.Kub*P.Mesh.Excite.Sub; zeros(NDofInt,P.Mesh.NExcite)];
% problem.Cu = [P.Model.Excite.Cub*P.Mesh.Excite.Sub; zeros(NDofInt,P.Mesh.NExcite)];
% problem.Mu = [P.Model.Excite.Mub*P.Mesh.Excite.Sub; zeros(NDofInt,P.Mesh.NExcite)];
% problem.M  = 0*blkdiag(P.Model.Rotor.M,zeros(NDofInt));
% problem.C  = 0*blkdiag(P.Model.Rotor.C,zeros(NDofInt));
% problem.K  = 0*blkdiag(P.Model.Rotor.K,zeros(NDofInt));
% problem.F0 = 0*[-P.Model.Fg + P.Model.Rotor.F0 - P.Model.Rotor.K*P.Model.x0; zeros(NDofInt,1)];

function Red = get_internal_indices(B,bCompressREB)
iDofIn = 0;
iInt = [];
iSum = {};
for i = 1:length(B)
    for j = 1:2
        if strcmp(B{i}.Model{j},'REB')
            REB = B{i}.Params{j};
            if bCompressREB
                iInt = [iInt; iDofIn+((1:REB.Model.NDof)-1)*REB.Setup.Z+1];
                for k = 1:REB.Model.NDof
                    iSum = [iSum; iDofIn+ (k-1)*REB.Setup.Z + (1:REB.Setup.Z)];
                end
                iDofIn = iDofIn + B{i}.NDofInt(j);
                Sparsity{i,j} = ones(B{i}.Params{j}.Model.NDof);
            else
                iInt = [iInt; iDofIn+(1:B{i}.NDofInt(j))];
                iSum = [iSum; num2cell(iDofIn+(1:B{i}.NDofInt(j)))];
                iDofIn  = iDofIn  + B{i}.NDofInt(j);
                Sparsity{i,j} = repmat(eye(REB.Setup.Z),B{i}.Params{j}.Model.NDof);
            end
        else
            %othwerise need all the internal states
            iInt = [iInt; iDofIn+(1:B{i}.NDofInt(j))];
            iSum = [iSum; num2cell(iDofIn+(1:B{i}.NDofInt(j)))];
            iDofIn  = iDofIn  + B{i}.NDofInt(j);
            Sparsity{i,j} = ones(B{i}.NDofInt(j));
        end
    end
end
Red.NDofInt = length(iInt);
Red.iInt = iInt;
Red.iSum = iSum;

Sparsity = Sparsity';
Red.Sparsity = blkdiag(Sparsity{:});
% Red.Sparsity = Red.Sparsity*0 + 1;