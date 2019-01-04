function x = rotor_equib(P,x0,O,A)
if nargin < 3 || isempty(O)
    O = 0;
end
if nargin < 4 || isempty(A)
    A = linspace(0,2*pi,100);
end
N = length(A);

if nargin < 2 || isempty(x0)
    x0 = [rand(P.Model.NDof,1)*1E-2;
          rand(P.Model.NDofInt*N,1)*1E-6];
else
    x0 = [x0(1:P.Model.NDof);
           repmat(x0(P.Model.NDof+1:end),N,1)];
end

Jb = get_sparsity(P.Bearing);
Jstr = [ones(P.Model.NDof)    ones(P.Model.NDof,P.Model.NDofInt*N);
        ones(P.Model.NDofInt*N,P.Model.NDof) kron(eye(N),Jb)];

options.maxit = 10;
options.ftol = 1E-8;
options.xtol = 1E-12;
options.jacobian = @rotor_equilibrium_jacob;
options.jacobianstructure = Jstr;
options.print_level = 0;

info.status = 2;

iter = 0;
bSuccess = 0;
while ~bSuccess && iter < 20
    [x,info] = fipopt([],x0,@rotor_equilibrium,options,P,O,A);
    constr = rotor_equilibrium(x,P,O,A);
    bSuccess = any(info.status == [0 1]);
    if ~bSuccess
        x0 = [rand(P.Model.NDof,1)*1E-2;
              rand(P.Model.NDofInt*N,1)*1E-6];
    else
        x0 = x;
    end
    iter = iter + 1;
end

J = full(rotor_equilibrium_jacob(x,P,O,A));
[xCG,xInt] = unpack_vector(x,P,N);
x = [xCG; 
    mean(xInt,2)];

if ~any(info.status == [0 1])
    error('Failed to find equilibirum position')
end

function [xCG,xInt] = unpack_vector(x,P,N)
xCG = x(1:P.Model.NDof);
xInt = reshape(x(P.Model.NDof+1:end),P.Model.NDofInt,N);

function constr = rotor_equilibrium(x,P,O,A)
wons = (0*A+1);
N = length(A);

[xCG,xInt] = unpack_vector(x,P,N);
uGnd = P.Mesh.Bearing.u0;     
     
States.O = O*wons; 
States.A = A;
States.x = (P.Model.A*xCG)*wons;
States.u = uGnd*wons;
States.xInt = xInt;
States.bSolve = 0;

Forces = bearingforces(P,States);

Fr  = P.Model.Rotor.K*xCG;
Fg  = P.Model.Fg;
Fb  = P.Model.A'*Forces.F;

constr = [Fg - Fr - mean(Fb,2); Forces.FInt(:)];

if any(isnan(constr))
    1
end

function J = rotor_equilibrium_jacob(x,P,O,A)
if 1
    N = length(A);
    wons = (0*A+1);

    [xCG,xInt] = unpack_vector(x,P,N);
    uGnd = P.Mesh.Bearing.u0;     

    States.O = O*wons; 
    States.A = A;               
    States.x = (P.Model.A*xCG)*wons;
    States.u = uGnd*wons;
    States.xInt = xInt;
    States.bSolve = 0;

    [Forces, Stiffness] = bearingforces(P,States);

    Fr  = P.Model.Rotor.K*xCG*wons;
    Fg  = P.Model.Fg*wons;
    Fb  = P.Model.A'*Forces.F;

    Jqq =  mtimesx(P.Model.A',mtimesx(Stiffness.Kqq,P.Model.A));
    Jqx =  mtimesx(P.Model.A',Stiffness.Kqx);
    Jxq =  mtimesx(Stiffness.Kxq,P.Model.A);
    Jxx =  Stiffness.Kxx;

    J = [-P.Model.Rotor.K-mean(Jqq,3) -catmat(Jqx,2)/N;
               catmat(Jxq,1)   blkmat(Jxx)];

else
    J = jacobian(@rotor_equilibrium,x,P,O,A);
end

if any(isnan(J))
    1
end

J = sparse(J);

function Sparsity = get_sparsity(B)
Sparsity = cell(length(B),2);
for i = 1:length(B)
    for j = 1:2
        if strcmp(B{i}.Model{j},'REB')
            REB = B{i}.Params{j};
            Sparsity{i,j} = repmat(eye(REB.Setup.Z),B{i}.Params{j}.Model.NDof);
        else
            Sparsity{i,j} = ones(B{i}.NDofInt(j));
        end
    end
end

Sparsity = Sparsity';
Sparsity = blkdiag(Sparsity{:});