function [V,d,W,S,R,O] = rotor_eig_cont(P,O,opt)
%ROTOT_EIG 

fun = @(x)statespace(x,P);
[V,d,W,O] = eigenshuffle_cont(fun,O,opt);

NDof = 2*size(P.A,2);
NSweep = length(O);
S = zeros(NDof,NDof,NSweep);
R = zeros(NDof,NDof,NSweep);

for i = 1:NSweep
    [R(:,:,i),S(:,:,i)] = statespace(O(i),P);
end

function [R,S] = statespace(O,P)
% Assemble it all into state-space form

    Z = 0*P.M;
    
    %Supply the matrices in the fixed domain
    S = [ P.M  O*P.G; %symmetric when Omega = 0
          Z     -P.K];
      
    R = -[P.C    P.K; %scaled by K so symmetric
          P.K     Z];
