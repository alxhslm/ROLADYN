function [Qi,Qo,Xr,Xz,wi,wo,K] = dynamic_contactlaw_harris(Contact,Geom,Race,Options,Fc,Fi,Fo,Ar,Az)

[wi,wo,Xr,Xz] = ball_newton(Contact,Geom,Race,Options,Fc,Fi,Fo,Ar,Az);
[dbi,dbo] = ball2contact(Geom,Ar,Az,Xr,Xz);

[Qi,Ki] = hertz_contactlaw(Contact.Inner.K,Contact.n,dbi,Contact.tol);
[Qo,Ko] = hertz_contactlaw(Contact.Outer.K,Contact.n,dbo,Contact.tol);

% if Options.bRaceCompliancei
%     [~,Kri] = race_compliance_loads(Race.Inner,-wi);
% else
    Kri = Inf;
% end
% if Options.bRaceComplianceo
%     [~,Kro] = race_compliance_loads(Race.Outer,wo);
% else
    Kro = Inf;
% end
K = 1./(1./Ki + 1./Ko + 1./Kri + 1./Kro);

function [dbi,dbo,ai,ao] = ball2contact(Geom,Ar,Az,Xr,Xz)
[Ai,Ao,ai,ao] = race_geometry(Xz,Xr,Az,Ar);
dbi = Ai - (Geom.RRacei-Geom.D/2);
dbo = Ao - (Geom.RRaceo-Geom.D/2);

function [wi,wo,Xr,Xz,iter] = ball_newton(Contact,Geom,Race,Options,Fc,Fi,Fo,Ar,Az)
A = sqrt(Az.^2 + Ar.^2);
Xz = (Az./A) .* ((Geom.RRaceo-Geom.D/2) + max(A - Geom.A0,0)/(1 + Contact.lambda));
Xr = (Ar./A) .* ((Geom.RRaceo-Geom.D/2) + max(A - Geom.A0,0)/(1 + Contact.lambda));

sz = size(A);
Xz = permute(Xz(:),[2 3 1]);
Xr = permute(Xr(:),[2 3 1]);
Az = permute(Az(:),[2 3 1]);
Ar = permute(Ar(:),[2 3 1]);
Fc = permute(Fc(:),[2 3 1]);
Fi = permute(Fi(:),[2 3 1]);
Fo = permute(Fo(:),[2 3 1]);
x = [Xr; Xz];

iter = 0;

% iFlex = [];
% if Options.bRaceCompliancei
%     iFlex(end+1) = 1;
% end
% if Options.bRaceComplianceo
%     iFlex(end+1) = 2;
% end
% if Options.bCentrifugal
%     iFlex(end+1) = 3;
% end
iFlex = [1;2];

% Qri = dn + NaN; Kri = dn + NaN;
% Qro = dn + NaN; Kro = dn + NaN;

dQ = x+Inf;

while any(abs(dQ(:))>1E-1)
    Xr = x(1,:,:);
    Xz = x(2,:,:);
    [dbi,dbo,ai,ao] = ball2contact(Geom,Ar,Az,Xr,Xz);
    
    [Qi,Ki] = hertz_contactlaw(Contact.Inner.K,Contact.n,dbi,Contact.tol);
    [Qo,Ko] = hertz_contactlaw(Contact.Outer.K,Contact.n,dbo,Contact.tol);

    [dbi_dXr,dbi_dXz] = ball_deriv(Ar-Xr,Az-Xz);
    dbi_dXr = -dbi_dXr; dbi_dXz = -dbi_dXz;
    [dai_dXr, dai_dXz] = angle_deriv(Ar-Xr,Az-Xz);
    dai_dXr = -dai_dXr; dai_dXz = -dai_dXz;
        
    [dao_dXr, dao_dXz] = angle_deriv(Xr,Xz);
    [dbo_dXr,dbo_dXz] = ball_deriv(Xr,Xz);
    
    Kdb = [Ki.*dbi_dXr.*sin(ai)-Ko.*dbo_dXr.*sin(ao) Ki.*dbi_dXz.*sin(ai)-Ko.*dbo_dXz.*sin(ao);
           Ki.*dbi_dXr.*cos(ai)-Ko.*dbo_dXr.*cos(ao) Ki.*dbi_dXz.*cos(ai)-Ko.*dbo_dXz.*cos(ao)];
        
    Kai = [(Qi.*cos(ai)+Fi.*sin(ai)).*dai_dXr  (Qi.*cos(ai)+Fi.*sin(ai)).*dai_dXz;
          (-Qi.*sin(ai)+Fi.*cos(ai)).*dai_dXr (-Qi.*sin(ai)+Fi.*cos(ai)).*dai_dXz];

    Kao = [(-Qo.*cos(ao)-Fo.*sin(ao)).*dao_dXr (-Qo.*cos(ao)-Fo.*sin(ao)).*dao_dXz;
            (Qo.*sin(ao)-Fo.*cos(ao)).*dao_dXr  (Qo.*sin(ao)-Fo.*cos(ao)).*dao_dXz];
    
    Kmat = Kdb + Kai + Kao;
    
    dQ = [Qi.*sin(ai) - Fi.*cos(ai) - Qo.*sin(ao) + Fo.*cos(ao);
          Qi.*cos(ai) + Fi.*sin(ai) - Qo.*cos(ao) - Fo.*sin(ao) + Fc];
    
    x(iFlex,:,:) = x(iFlex,:,:) - mtimesx(minvx(Kmat(iFlex,iFlex,:)),dQ(iFlex,:,:));
    
    iter = iter + 1;
    if iter > 100
        break
    end
end

Xr = reshape(x(1,:,:),sz);
Xz = reshape(x(2,:,:),sz);

wi = 0*Xr;
wo = 0*Xz;

function [dbo_dXr,dbo_dXz] = ball_deriv(Xr,Xz)
A = sqrt(Xr.^2 + Xz.^2);
dbo_dXr = Xr./A;
dbo_dXz = Xz./A;

function [dao_dXr, dao_dXz] = angle_deriv(Xr,Xz)
y = Xz./Xr;
dao_dXr = 1./(1+y.^2) .* (-Xz./Xr.^2);
dao_dXz = 1./(1+y.^2) .* (1./Xr);