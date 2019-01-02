function [R,S,A] = rotor_ss(O,FE)
% Assemble it all into state-space form

if length(O) > 1
    for i = 1:length(O)
        [R(:,:,i),S(:,:,i),A(:,:,i)] = rotor_ss(O(i),FE);
    end
    return;
end

%Supply the matrices in the fixed domain
C = FE.C+O*FE.G;
M = FE.M;
K = FE.K;

Z = 0*M;
I = eye(size(M,1));

if rank(M)<size(M,1)
    B = null(M,'r');
    A = null(B'*K,'r');
else
    A = I;
end

bGenerlised = 1;

% N = null(FE.M);
% Q = orth(FE.M);
% 
% Mmm = Q'*FE.M*Q;
% Cmm = Q'*Ctot*Q;
% Kmm = Q'*FE.K*Q;
% 
% Cms = Q'*Ctot*N;
% Kms = Q'*FE.K*N;
% 
% Csm = N'*Ctot*Q;
% Ksm = N'*FE.K*Q;
% 
% Css = N'*Ctot*N;
% Kss = N'*FE.K*N;
% 
% Zmm = 0*Kmm;
% Zms = 0*Kms;
% Zsm = 0*Ksm;
% Imm = eye(size(Zmm));
% 
% S0 = [Mmm    Zmm  Cms;
%       Zmm   -Imm  Zms;
%       Zsm    Zsm  Css];
%  
% R0 = -[Cmm  Kmm  Kms;
%        Zmm  Imm  Zms;
%        Csm  Ksm  Kss];
% 
% A = [ Q  0*Q 0*N;
%      0*Q  Q   N];
% 
% S2 = A'*S*A;
% R2 = A'*R*A;

if bGenerlised
    %Generalised eigenvalue problem
    S = [  M    C; %symmetric when Omega = 0
           Z     I];

    R = -[Z   K; %scaled by K so symmetric
         -I   Z];
else
    %Standard eigenvalue problem
    S = [ I  Z; %symmetric when Omega = 0
          Z  I];

    R = [-M\C   -M\K; %scaled by K so symmetric
           I      Z];
end

% A = eye(size(S,1));
% if rank(S) < size(S,2)
%     B = null(S);
%     A = null(B'*R);
% end
% 
% S2 = A'*S*A;
% R2 = A'*R*A;
% if rank(S2) < size(S2,2)
%     B = null(S2);
%     A = A*null(B'*R2);
% end
