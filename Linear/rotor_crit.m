function [V,d,W,S,R]  = rotor_crit(FE,N)

if nargin < 2
    N = 1;
end

NDof = 2*size(FE.A,2);
NSweep = length(N);

S = zeros(NDof,NDof,NSweep);
R = zeros(NDof,NDof,NSweep);

for i = 1:NSweep
    [R(:,:,i),S(:,:,i)] = statespace(N(i),FE);
end
[V,d,W] = eigenshuffle(R,S,N);

function [R,S] = statespace(n,FE)
Z = 0*FE.M;

S = [(FE.M*n^2-1i*n*FE.G)   Z; %symmetric when Omega = 0
          Z      -FE.K];
      
R = -[n*FE.C   FE.K; %scaled by K so symmetric
      FE.K    Z];