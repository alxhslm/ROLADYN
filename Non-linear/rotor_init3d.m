function [X,x,t] = rotor_init3d(hbm,P,w,A)
if hbm.options.bUseStandardHBM
    [X,x,t] = rotor_init(hbm,P,w,A);
    return;
end
w0 = w*hbm.harm.rFreqRatio + hbm.harm.wFreq0;
O = w0(1);
wc = w0(2);

NDof = P.Model.NDof;
NHarm = hbm.harm.NHarm;
NFreq = hbm.harm.NFreq;
Nfft = hbm.harm.Nfft;

U = A*rotor_hbm_excite(hbm,problem,w0)
  
%work out frequencies
w = hbm.harm.kHarm(:,1)*w0(1) + hbm.harm.kHarm(:,2)*w0(2);

M = P.Model.M;
C = P.Model.C;
K = P.Model.K;
Mu = P.Model.Excite.Me;
Cu = P.Model.Excite.Ce;
Ku = P.Model.Excite.Ke;

% [M,C,K,Mu,Cu,Ku] = getMatrices(hbm,problem,[P.Model.x0;P.Model.xInt],w0);

clear Fe
Fe(1,1:NDof) = P.Model.Fg;
for i = 2:NFreq
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

switch hbm.options.aft_method
    case 'fft'
        %create the time series from the fourier series
        xCG     = freq2time3d(X,NHarm,hbm.harm.iSub,Nfft)';
        xdotCG  = freq2time3d(Xdot,NHarm,hbm.harm.iSub,Nfft)';
        xddotCG = freq2time3d(Xddot,NHarm,hbm.harm.iSub,Nfft)';
        
        fe = freq2time3d(Fe,NHarm,hbm.harm.iSub,Nfft)';
    case 'mat'
        %create the time series from the fourier series
        xCG     = real(hbm.nonlin.IFFT*X)';
        xdotCG  = real(hbm.nonlin.IFFT*Xdot)';
        xddotCG = real(hbm.nonlin.IFFT*Xddot)';

        fe = real(hbm.nonlin.IFFT*Fe)';
end

t1  = (0:Nfft(1)-1)/Nfft(1)*2*pi/w0(1);
t2 = (0:Nfft(2)-1)/Nfft(2)*2*pi/w0(2);
[t1,t2] = ndgrid(t1,t2);
t = [t1(:)'; t2(:)'];

%get the internal states
x0 = P.Model.x0*(0*t(1,:)+1);
u0 = P.Mesh.Bearing.u0*(0*t(1,:)+1);

%compute the response of the linearised system
Fgyro = O*P.Model.Rotor.G*xdotCG;
fR = P.Model.Rotor.K*(xCG-x0)   + P.Model.Rotor.C*xdotCG   + P.Model.Rotor.M*xddotCG      + P.Model.Rotor.F0;
fB = P.Model.Bearing.K*(xCG-x0) + P.Model.Bearing.C*xdotCG + P.Model.Bearing.F0 + Fgyro;

%now compute the response of the fully non-linear system
States.x     = P.Model.A*xCG;
States.xdot  = P.Model.A*xdotCG;
States.xddot = P.Model.A*xddotCG;
States.u     = uGnd;
States.udot  = udotGnd;
States.uddot = uddotGnd;
States.A = O*t(P.Model.iRot,:);
States.O = O + t(P.Model.iRot,:)*0;
States.bSolve = 1;
Forces = bearingforces(P,States);
fNL = P.Model.A'*Forces.F + Fgyro;
xInt = Forces.xInt;

%comparison of rotor states
iPlot = [1 2];
NPlot = length(iPlot);
figure
if Nfft(2) > 1
    for i = 1:NPlot
        subplot(NPlot,1,i)
        surf(t2,t1,reshape(xCG(iPlot(i),:),Nfft(1),Nfft(2)))
    end
else
    plot(t1,xCG(iPlot,:))
end
xlabel('t_2')
ylabel('t_1')

%comparison of bearing forces
figure
if Nfft(2) > 1
    for i = 1:NPlot
        subplot(NPlot,1,i)
        surf(t2,t1,reshape(fNL(iPlot(i),:),Nfft(1),Nfft(2)),'EdgeColor','none')
        hold on
        mesh(t2,t1,reshape(fB(iPlot(i),:),Nfft(1),Nfft(2)),'FaceColor','none')
    end
    xlabel('t_2')
    ylabel('t_1')
    zlabel('f')
else
    for i = 1:NPlot
        subplot(NPlot,1,i)
        plot(t1,fNL(iPlot(i),:));
        hold on
        plot(t1,fB(iPlot(i),:),'o');
    end
    xlabel('t_1')
    ylabel('f')
end

if ~isempty(xInt)
    figure
    hold on
    if Nfft(2) > 1
        subplot(2,1,1)
        hold on
        for i = 1:size(xInt,1)/2
            surf(t2,t1,reshape(abs(xInt(i,:)),Nfft(1),Nfft(2)))
        end
        subplot(2,1,2)
        hold on
        for i = (size(xInt,1)/2+1):size(xInt,1)
            surf(t2,t1,reshape(xInt(i,:),Nfft(1),Nfft(2)))
        end

        xlabel('t_2')
        ylabel('t_1')
        zlabel('xInt')
    else
        for i = 1:size(xInt,1)
            plot(t1,xInt(i,:));
        end
        xlabel('t_1')
        ylabel('f')
    end
end

% FNL = time2freq3d(fNL',NHarm,hbm.harm.iSub,Nfft);   
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

XInt = time2freq3d(xInt(P.Model.Reduced.iInt,:)',NHarm,hbm.harm.iSub,Nfft);    

% XInt2 = hbm.nonlin.FFT*xInt';
x = [xCG; xInt];
X = [X XInt];

if Nfft(2) > 1
    x = reshape(x',Nfft(1),Nfft(2),[]);
    t = reshape(t',Nfft(1),Nfft(2),[]);
end

X = packdof(X,hbm.harm);