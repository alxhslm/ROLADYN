function P = rotor_equib(P,O,A)
if nargin < 2 || isempty(O)
    O = 0;
end
if nargin < 3 || isempty(A)
    A = linspace(0,2*pi,100);
end
N = length(A);

if isfield(P.Mesh,'x0') && length(P.Model.x0) == P.Model.NDof
    x0 = [P.Mesh.x0;
          repmat(P.Mesh.xInt,N,1)];
else
    x0 = [rand(P.Model.NDof,1)*1E-2;
          rand(P.Model.NDofInt*N,1)*1E-6];
end

Jb = get_sparsity(P.Bearing);
Jstr = [ones(P.Model.NDof)    ones(P.Model.NDof,P.Model.NDofInt*N);
        ones(P.Model.NDofInt*N,P.Model.NDof) kron(eye(N),Jb)];

options.maxit = 10;
options.ftol = 1E-8;
options.xtol = 1E-12;
options.jacobian = @rotor_equilibrium_jacob;
options.jacobianstructure = Jstr;
options.print_level = 5;

info.status = 2;

iter = 0;
bSuccess = 0;
while ~bSuccess && iter < 20
    [x,info] = fipopt([],x0,@rotor_equilibrium,options,P,O,A);
    bSuccess = any(info.status == [0 1]);
    x0 = x;
    iter = iter + 1;
end

[xCG,xInt] = unpack_vector(x,P,N);

P.Model.x0   = xCG;
P.Model.xInt = mean(xInt,2);

P.Mesh.x0   = P.Model.A * P.Model.x0;
P.Mesh.xInt = P.Model.xInt;

wons = (0*A+1);
States.O = O*wons; 
States.A = A;
States.x = (P.Model.Bearing.S*xCG)*wons;
States.xInt = xInt;
States.bSolve = 0;

[Forces, Stiffness] = bearingforces(P,States);

P.Model.Rotor.F0    = P.Model.Rotor.K*P.Model.x0;
P.Model.Stator.F0   = P.Model.Stator.K*P.Model.x0;
P.Model.Bearing.F0  = P.Model.Bearing.S'*mean(Forces.F,2);
P.Model.Bearing.K   = P.Model.Bearing.S'*mean(Stiffness.Kqq,3)*P.Model.Bearing.S;
P.Model.K           = P.Model.Rotor.K + P.Model.Stator.K + P.Model.Bearing.K;

P.Mesh.Rotor.F0    = P.Mesh.Rotor.K*P.Mesh.x0;
P.Mesh.Stator.F0   = P.Mesh.Stator.K*P.Mesh.x0;
P.Mesh.Bearing.F0  = P.Mesh.Bearing.S'*mean(Forces.F,2);
P.Mesh.Bearing.Fb  = mean(Forces.F,2);
P.Mesh.Bearing.K   = P.Mesh.Bearing.S'*mean(Stiffness.Kqq,3)*P.Mesh.Bearing.S;
P.Mesh.Bearing.Kb  = mean(Stiffness.Kqq,3);
P.Mesh.K           = P.Mesh.Rotor.K + P.Mesh.Stator.K + P.Mesh.Bearing.K;

if ~bSuccess
    error('Failed to find equilibrium position')
end

function [xCG,xInt] = unpack_vector(x,P,N)
xCG = x(1:P.Model.NDof);
xInt = reshape(x(P.Model.NDof+1:end),P.Model.NDofInt,N);

function constr = rotor_equilibrium(x,P,O,A)
wons = (0*A+1);
N = length(A);

[xCG,xInt] = unpack_vector(x,P,N);
     
States.O = O*wons; 
States.A = A;
States.x = (P.Model.Bearing.S*xCG)*wons;
States.xInt = xInt;
States.bSolve = 0;

Forces = bearingforces(P,States);

Fr  = (P.Model.Rotor.K+P.Model.Stator.K)*xCG;
Fg  = P.Model.Fg;
Fb  = P.Model.Bearing.S'*Forces.F;

constr = [Fg - Fr - mean(Fb,2); Forces.FInt(:)];

function J = rotor_equilibrium_jacob(x,P,O,A)
N = length(A);
wons = (0*A+1);

[xCG,xInt] = unpack_vector(x,P,N);

States.O = O*wons; 
States.A = A;               
States.x = (P.Model.Bearing.S*xCG)*wons;
States.xInt = xInt;
States.bSolve = 0;

[Forces, Stiffness] = bearingforces(P,States);

Jqq =  mtimesx(P.Model.Bearing.S',mtimesx(Stiffness.Kqq,P.Model.Bearing.S));
Jqx =  mtimesx(P.Model.Bearing.S',Stiffness.Kqx);
Jxq =  mtimesx(Stiffness.Kxq,P.Model.Bearing.S);
Jxx =  Stiffness.Kxx;

Klin = P.Model.Rotor.K+P.Model.Stator.K;
J = [-Klin-mean(Jqq,3) -catmat(Jqx,2)/N;
        catmat(Jxq,1)    blkmat(Jxx)];

J = sparse(J);

function Sparsity = get_sparsity(B)
Sparsity = cell(length(B),2);
for i = 1:length(B)
    if strcmp(B{i}.Model,'REB')
        REB = B{i}.Params;
        Sparsity{i} = repmat(eye(REB.Setup.Z),B{i}.Params.Model.NDof);
    else
        Sparsity{i} = ones(B{i}.NDofInt);
    end
end

Sparsity = Sparsity';
Sparsity = blkdiag(Sparsity{:});