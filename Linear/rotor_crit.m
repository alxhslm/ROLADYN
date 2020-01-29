function [V,d,W,S,R]  = rotor_crit(FE,N)

if nargin < 2
    N = 1;
end

NDof = 2*size(FE.A,2);
NSweep = length(N);

S = zeros(NDof,NDof,NSweep);
R = zeros(NDof,NDof,NSweep);

for i = 1:NSweep
    [R(:,:,i),S(:,:,i),A(:,:,i)] = statespace(N(i),FE);
    A2(:,:,i) = blkdiag(A(:,:,i),A(:,:,i));
end

At = mtransposex(A2);
S2 = mtimesx(At,mtimesx(S,A2));
R2 = mtimesx(At,mtimesx(R,A2));

[V,d,W] = eigenshuffle(R2,S2,N);

for i = 1:NSweep
    V2(:,:,i) = A2(:,:,i)*V(:,:,i);
    W2(:,:,i) = A2(:,:,i)*W(:,:,i);
end
V = V2;
W = W2;

function [R,S,A] = statespace(n,FE)
Z = 0*FE.M;
I = eye(size(FE.M,1));

if rank(FE.M)<size(FE.M,1)
    B = null(FE.M,'r');
    A = null(B'*FE.K,'r');
else
    A = I;
end


S = [(FE.M*n^2-1i*n*FE.G)   Z; %symmetric when Omega = 0
          Z      -FE.K];
      
R = -[n*FE.C   FE.K; %scaled by K so symmetric
      FE.K    Z];