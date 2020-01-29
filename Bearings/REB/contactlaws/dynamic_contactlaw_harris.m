function [Qi,Qo,Xr,Xz,wi,wo,Ktoti,Ktoto] = dynamic_contactlaw_harris(Contact,Geom,Race,Options,Fc,Fi,Fo,Ar,Az)

[wi,wo,Xr,Xz] = ball_newton(Contact,Geom,Race,Options,Fc,Fi,Fo,Ar,Az);
[Ai,Ao,ai,ao] = race_geometry(Xz,Xr,Az,Ar);
dbi = Ai - (Geom.RRacei-Geom.D/2) - wi;
dbo = Ao - (Geom.RRaceo-Geom.D/2) - wo;

[Qi,Ki] = hertz_contactlaw(Contact.Inner.K,Contact.n,dbi,Contact.tol);
[Qo,Ko] = hertz_contactlaw(Contact.Outer.K,Contact.n,dbo,Contact.tol);

if Options.bRaceCompliancei
    [~,Kri] = race_compliance_loads(Race.Inner,-wi);
else
    Kri = Inf;
end
if Options.bRaceComplianceo
    [~,Kro] = race_compliance_loads(Race.Outer,wo);
else
    Kro = Inf;
end
Ktoti = 1./(1./Ki + 1./Kri + 1./(Ko.*cos(ai-ao).^2) + 1./(Kro.*cos(ai-ao).^2));
Ktoto = 1./(1./Ko + 1./Kro + 1./(Ki.*cos(ai-ao).^2) + 1./(Kri.*cos(ai-ao).^2));

function [wi,wo,Xr,Xz,iter] = ball_newton(Contact,Geom,Race,Options,Fc,Fi,Fo,Ar,Az)
if Options.bRaceCompliancei
    Kri = Race.Inner.K;
else
    Kri = Inf;
end
if Options.bRaceCompliancei
    Kro = Race.Outer.K;
else
    Kro = Inf;
end

%sensible initial guesses
A = sqrt(Az.^2 + Ar.^2);
Xz = (Az./A) .* ((Geom.RRaceo-Geom.D/2) + max(A - Geom.A0,0)/(1 + Contact.lambda));
Xr = (Ar./A) .* ((Geom.RRaceo-Geom.D/2) + max(A - Geom.A0,0)/(1 + Contact.lambda));
Q0 = hertz_contactlaw(Contact.K,Contact.n,A - Geom.A0,Contact.tol);
wi = Q0/Kri;
wo = Q0/Kro;

%find elements out of contact
dn = A - Geom.A0;
dn_crit = (Fc/Contact.Outer.K).^(1/Contact.n);
iOuterOnly = dn<dn_crit;
wi(iOuterOnly) = 0;
wo(iOuterOnly) = (Fc(iOuterOnly)/Kro);
Xr(iOuterOnly) = (Geom.RRaceo-Geom.D/2) + dn_crit(iOuterOnly);
Xz(iOuterOnly) = 0;

%make everything 1D
sz = size(A);
Az = permute(Az(:),[2 3 1]);
Ar = permute(Ar(:),[2 3 1]);
Fc = permute(Fc(:),[2 3 1]);
Fi = permute(Fi(:),[2 3 1]);
Fo = permute(Fo(:),[2 3 1]);

wi = permute(wi(:),[2 3 1]);
wo = permute(wo(:),[2 3 1]);
Xr = permute(Xr(:),[2 3 1]);
Xz = permute(Xz(:),[2 3 1]);

iFlex = [];
if Options.bRaceCompliancei
    iFlex(end+1) = 1;
end
if Options.bRaceComplianceo
    iFlex(end+1) = 2;
end
iFlex = [iFlex 3 4];

x = [wi; wo; Xr; Xz];
dQ = x + Inf;
iter = 0;
Z = 0*Az;

f = @(y)myfunc(y,Contact,Geom,Race,Options,Ar,Az,Fc,Fi,Fo);
% dQ = feval(f,x);

while any(abs(dQ(:))>1E-3)
%     wi = x(1,:,:);
%     wo = x(2,:,:);
%     Xr = x(3,:,:);
%     Xz = x(4,:,:);
    
    dQ = feval(f,x);
    
    h = 1E-10;
    for i = 1:4
        x2 = x;
        x2(i,:) = x2(i,:) + h;
        Kmat(:,i,:) = permute(feval(f,x2)-dQ,[1 3 2])/h;
    end
    
%     [Ai,Ao,ai,ao] = race_geometry(Xz,Xr,Az,Ar);
%     dbi = Ai - (Geom.RRacei-Geom.D/2) - wi;
%     dbo = Ao - (Geom.RRaceo-Geom.D/2) - wo;
% 
%     [Qi,Ki] = hertz_contactlaw(Contact.Inner.K,Contact.n,dbi,Contact.tol);
%     [Qo,Ko] = hertz_contactlaw(Contact.Outer.K,Contact.n,dbo-wo,Contact.tol);
% %     [Q,K] = hertz_contactlaw(Contact.K,Contact.n,A - Geom.A0 -wo-wi,Contact.tol);
%      
%     if Options.bRaceCompliancei
%         [Qri,Kri] = race_compliance_loads(Race.Inner,-wi); Qri = -Qri;
%     else
%         Qri = Qi;
%         Kri = Az+Inf;
%     end
%     if Options.bRaceComplianceo
%         [Qro,Kro] = race_compliance_loads(Race.Outer,wo);
%     else
%         Qro = Qo;
%         Kro = Az+Inf;
%     end
%       
%     dbi_dXr = -cos(ai);
%     dbi_dXz = -sin(ai);
%     
%     dai_dXr =  sin(ao)./Ai;
%     dai_dXz = -cos(ao)./Ai;
% 
%     dbo_dXr = cos(ao);
%     dbo_dXz = sin(ao);
%     
%     dao_dXr = -sin(ao)./Ao;
%     dao_dXz =  cos(ao)./Ao;
%     
%     Kdbi = [Ki.*dbi_dXr.*sin(ai) Ki.*dbi_dXz.*sin(ai);
%             Ki.*dbi_dXr.*cos(ai) Ki.*dbi_dXz.*cos(ai)];
%        
%     Kdbo =-[Ko.*dbo_dXr.*sin(ao) Ko.*dbo_dXz.*sin(ao);
%             Ko.*dbo_dXr.*cos(ao) Ko.*dbo_dXz.*cos(ao)];
% 
%     Kai = [(Qi.*cos(ai)+Fi.*sin(ai)).*dai_dXr  (Qi.*cos(ai)+Fi.*sin(ai)).*dai_dXz;
%           (-Qi.*sin(ai)+Fi.*cos(ai)).*dai_dXr (-Qi.*sin(ai)+Fi.*cos(ai)).*dai_dXz];
% 
%     Kao = [(-Qo.*cos(ao)-Fo.*sin(ao)).*dao_dXr (-Qo.*cos(ao)-Fo.*sin(ao)).*dao_dXz;
%             (Qo.*sin(ao)-Fo.*cos(ao)).*dao_dXr  (Qo.*sin(ao)-Fo.*cos(ao)).*dao_dXz];
% 
%     Krace =-[Ki+Kri   Z;    
%                 Z   Ko+Kro];
%     Kcent = Kdbi + Kai + Kdbo + Kao;
%     Kcross = [Ki.*dbi_dXr Ki.*dbi_dXz;
%               Ko.*dbo_dXr Ko.*dbo_dXz];
%           
%     Kcross2 = -[Ki.*sin(ai) -Ko.*sin(ao);
%                 Ki.*cos(ai) -Ko.*cos(ao)];
%           
%     Kmat = [  Krace      Kcross;
%               Kcross2    Kcent];
% 
%     dQ = [Qi-Qri;
%           Qo-Qro;
%           Qi.*sin(ai) - Fi.*cos(ai) - Qo.*sin(ao) + Fo.*cos(ao);
%           Qi.*cos(ai) + Fi.*sin(ai) - Qo.*cos(ao) - Fo.*sin(ao) + Fc];
    
    x(iFlex,:,:) = x(iFlex,:,:) - mtimesx(minvx(Kmat(iFlex,iFlex,:)),dQ(iFlex,:,:));
    
    iter = iter + 1;
    if iter > 100
        break
    end
end

wi = reshape(x(1,:,:),sz);
wo = reshape(x(2,:,:),sz);
Xr = reshape(x(3,:,:),sz);
Xz = reshape(x(4,:,:),sz);

% function [dbo_dXr,dbo_dXz] = ball_deriv(ai)
% dbo_dXr = cos(ai);
% dbo_dXz = sin(ai);
% 
% function [dai_dXr, dai_dXz] = angle_deriv(ai,X)%Xr,Xz)
% % y = Xz./Xr;
% dai_dXr = -sin(ai)./X;%1./(1+y.^2) .* (-Xz./Xr.^2);
% dai_dXz =  cos(ai)./X;%1./(1+y.^2) .* (1./Xr);

function dQ = myfunc(x,Contact,Geom,Race,Options,Ar,Az,Fc,Fi,Fo)
wi = x(1,:,:);
wo = x(2,:,:);
Xr = x(3,:,:);
Xz = x(4,:,:);

[Ai,Ao,ai,ao] = race_geometry(Xz,Xr,Az,Ar);
dbi = Ai - (Geom.RRacei-Geom.D/2) - wi;
dbo = Ao - (Geom.RRaceo-Geom.D/2) - wo;

Qi = hertz_contactlaw(Contact.Inner.K,Contact.n,dbi,Contact.tol);
Qo = hertz_contactlaw(Contact.Outer.K,Contact.n,dbo,Contact.tol);

if Options.bRaceCompliancei
    Qri = race_compliance_loads(Race.Inner,-wi); Qri = -Qri;
else
    Qri = Qi;
end
if Options.bRaceComplianceo
    Qro = race_compliance_loads(Race.Outer,wo);
else
    Qro = Qo;
end

dQ = [Qi-Qri;
      Qo-Qro;
      Qi.*sin(ai) - Fi.*cos(ai) - Qo.*sin(ao) + Fo.*cos(ao);
      Qi.*cos(ai) + Fi.*sin(ai) - Qo.*cos(ao) - Fo.*sin(ao) + Fc];