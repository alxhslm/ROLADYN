function [V,d,W,S,R] = rotor_eig_async(P,O1,O2)
%ROTOT_EIG 

NDof = 2*size(P.A,2);

S = zeros(NDof,NDof,length(O1),length(O2));
R = zeros(NDof,NDof,length(O1),length(O2));

for i = 1:length(O1)
    for j = 1:length(O2)
        [R(:,:,i,j),S(:,:,i,j)] = statespace(O1(i),O2(j),P);
    end
end
[V,d,W] = eigenshuffle_async(R,S,O1,O2);

function [R,S] = statespace(O1,O2,P)
% Assemble it all into state-space form

    Z = 0*P.M;
    G = P.G;
    N = size(P.M,1);
    G(1:N/2,1:N/2) = G(1:N/2,1:N/2) * O1;
    G(N/2 + (1:N/2),N/2 + (1:N/2)) = G(N/2 + (1:N/2),N/2 + (1:N/2)) * O2;
    
    %Supply the matrices in the fixed domain
    S = [ P.M    G; %symmetric when Omega = 0
          Z     -P.K];
      
    R = -[P.C    P.K; %scaled by K so symmetric
          P.K     Z];
