function [X,x,t] = rotor_init(hbm,P,w0,A)
O = w0;

NDof  = P.Model.NDof;
NHarm = hbm.harm.NHarm(1);
NFreq = hbm.harm.NFreq;
Nfft  = hbm.harm.Nfft(1);

Fe = excitation_hbm(P);
U = A*packharm(Fe{1},NHarm);

%work out frequencies
w = (0:NHarm)'*w0;

M = P.Model.M;
C = P.Model.C;
K = P.Model.K;
Mu = P.Model.Excite.Me;
Cu = P.Model.Excite.Ce;
Ku = P.Model.Excite.Ke;

% [M,C,K,Mu,Cu,Ku] = getMatrices(hbm,problem,[P.Model.x0;P.Model.xInt],w0);

clear Fe
for i = 1:NFreq
    Fe(i,:) = ((Ku + 1i*w(i)*Cu - w(i)^2*Mu)*U(i,:).').';
end

X(1,1:NDof) = P.Model.x0;
for i = 2:NFreq
    X(i,:) = ((K + 1i*w(i)*C - w(i)^2*M)\Fe(i,:).').';
end

%compute the fourier coefficients of the derivatives
Wx = repmat(1i*w,1,size(X,2));
Xdot  = X.*Wx;
Xddot = Xdot.*Wx;

%precompute the external inputs
Wu = repmat(1i*w,1,size(U,2));
Udot  = U.*Wu;
Uddot = Udot.*Wu;

%create the time series from the fourier series
xCG     = freq2time(X,NHarm,Nfft)';
xdotCG  = freq2time(Xdot,NHarm,Nfft)';
xddotCG = freq2time(Xddot,NHarm,Nfft)';

%create the vector of inputs
uGnd     = freq2time(U*P.Mesh.Excite.Sgd.',NHarm,Nfft)';
udotGnd  = freq2time(Udot*P.Mesh.Excite.Sgd.',NHarm,Nfft)';
uddotGnd = freq2time(Uddot*P.Mesh.Excite.Sgd.',NHarm,Nfft)';

fe = freq2time(Fe,NHarm,Nfft)';

Fgyro = O*P.Model.Rotor.G*xdotCG;

t = (0:Nfft-1)/Nfft*2*pi/w0;
A = O*t;

O = O + 0*A;

%get the internal states
x0 = P.Model.x0*(0*O+1);
u0 = P.Mesh.Bearing.u0*(0*O+1);

fe = fe + P.Model.Fg;
fL = P.Model.Rotor.K*(xCG-x0)   + P.Model.Rotor.C*xdotCG   + P.Model.Rotor.M*xddotCG + P.Model.Rotor.F0;
fB = P.Model.Bearing.K*(xCG-x0) + P.Model.Bearing.C*xdotCG + P.Model.Excite.Kgd*(uGnd-u0) + P.Model.Bearing.F0 + Fgyro;

States.x     = P.Model.A*xCG;
States.xdot  = P.Model.A*xdotCG;
States.xddot = P.Model.A*xddotCG;
States.u     = uGnd;
States.udot  = udotGnd;
States.uddot = uddotGnd;
States.A = A;
States.O = O;
States.bSolve = 1;
Forces = bearingforces(P,States);
fNL = P.Model.A'*Forces.F + Fgyro;
xInt = Forces.xInt;


% FNL = time2freq(fNL',NHarm,Nfft);   
% FNL2 = hbm.nonlin.FFT*fNL'; 
% Nf = size(fNL,1);
% 
% figure
% for i = 1:Nf
%     subplot(Nf,1,i)   
%     plot(-fNL(i,:));
%     hold on
%     plot(-fB(i,:));
% end
% 
% figure
% for i = 1:Nf
%     subplot(Nf,1,i)   
%     plot(fL(i,:)-fe(i,:));
%     hold on
%     plot(-fB(i,:));
% end

XInt = time2freq(xInt',NHarm,Nfft);

% XInt2 = hbm.aft.FFT*xInt';
x = [xCG; xInt];
X = [X XInt];