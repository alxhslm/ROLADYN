function [Qi,Qo,Xr,Xz,wi,wo] = dynamic_contactlaw_harris(B,Oi,Oo,Ar,Az)

[Qi,Qo,wi,wo,Xr,Xz] = ball_newton(B,Oi,Oo,Ar,Az);

function [Qi,Qo,wi,wo,Xr,Xz,iter] = ball_newton(B,Oi,Oo,Ar,Az)
if B.Options.bRaceCompliancei
    Kri = Race.Inner.K;
else
    Kri = Inf;
end
if B.Options.bRaceCompliancei
    Kro = Race.Outer.K;
else
    Kro = Inf;
end

%sensible initial guesses
A = sqrt(Az.^2 + Ar.^2);
Xz = (Az./A) .* ((B.Geometry.RRaceo-B.Geometry.D/2) + max(A - B.Geometry.A0,0)/(1 + B.Contact.lambda));
Xr = (Ar./A) .* ((B.Geometry.RRaceo-B.Geometry.D/2) + max(A - B.Geometry.A0,0)/(1 + B.Contact.lambda));
Q0 = hertz_contactlaw(B.Contact.K,B.Contact.n,A - B.Geometry.A0,B.Contact.tol);
wi = Q0/Kri;
wo = Q0/Kro;

%find elements out of B.Contact
dn = A - B.Geometry.A0;
Fc = dynamic_ball_loads(B,0*A,0*A,Oi,Oo);
dn_crit = (Fc/B.Contact.Outer.K).^(1/B.Contact.n);
iOuterOnly = dn<dn_crit;
wi(iOuterOnly) = 0;
wo(iOuterOnly) = (Fc(iOuterOnly)/Kro);
Xr(iOuterOnly) = (B.Geometry.RRaceo-B.Geometry.D/2) + dn_crit(iOuterOnly);
Xz(iOuterOnly) = 0;

%make everything 1D
sz = size(Ar);

Az = permute(Az(:),[2 3 1]);
Ar = permute(Ar(:),[2 3 1]);

Oi = permute(Oi(:),[2 3 1]);
Oo = permute(Oo(:),[2 3 1]);

wi = permute(wi(:),[2 3 1]);
wo = permute(wo(:),[2 3 1]);
Xr = permute(Xr(:),[2 3 1]);
Xz = permute(Xz(:),[2 3 1]);

iFlex = [];
if B.Options.bRaceCompliancei
    iFlex(end+1) = 1;
end
if B.Options.bRaceComplianceo
    iFlex(end+1) = 2;
end
iFlex = [iFlex 3 4];

x = [wi; wo; Xr; Xz];
dQ = x + Inf;
iter = 0;

while any(abs(dQ(:))>1E-3)   
    [J,dQ] = ball_equib_jacob(x,B,Oi,Oo,Ar,Az);
    x(iFlex,:,:) = x(iFlex,:,:) - mtimesx(minvx(J(iFlex,iFlex,:)),dQ(iFlex,:,:));
    
    iter = iter + 1;
    if iter > 100
        break
    end
end

[~,Qi,Qo] = ball_equib(x,B,Oi,Oo,Ar,Az);

wi = reshape(x(1,:,:),sz);
wo = reshape(x(2,:,:),sz);
Xr = reshape(x(3,:,:),sz);
Xz = reshape(x(4,:,:),sz);

Qi = reshape(Qi,sz);
Qo = reshape(Qo,sz);

function [Kmat,dQ] = ball_equib_jacob(x,B,Oi,Oo,Ar,Az)
dQ = ball_equib(x,B,Oi,Oo,Ar,Az);

h = 1E-10;
Kmat = zeros(size(dQ,1),size(x,1),size(x,3));
for i = 1:4
    x2 = x;
    x2(i,:) = x2(i,:) + h;
    dQ2 = ball_equib(x2,B,Oi,Oo,Ar,Az);
    Kmat(:,i,:) = permute(dQ2-dQ,[1 3 2])/h;
end

function [dQ,Qi,Qo] = ball_equib(x,B,Oi,Oo,Ar,Az)
wi = x(1,:,:);
wo = x(2,:,:);
Xr = x(3,:,:);
Xz = x(4,:,:);

[Ai,Ao,ai,ao] = race_geometry(Xz,Xr,Az,Ar);
dbi = Ai - (B.Geometry.RRacei-B.Geometry.D/2) - wi;
dbo = Ao - (B.Geometry.RRaceo-B.Geometry.D/2) - wo;

Qi = hertz_contactlaw(B.Contact.Inner.K,B.Contact.n,dbi,B.Contact.tol);
Qo = hertz_contactlaw(B.Contact.Outer.K,B.Contact.n,dbo,B.Contact.tol);

if B.Options.bRaceCompliancei
    Qri = race_compliance_loads(Race.Inner,-wi); Qri = -Qri;
else
    Qri = Qi;
end
if B.Options.bRaceComplianceo
    Qro = race_compliance_loads(Race.Outer,wo);
else
    Qro = Qo;
end

[Fc,Fi,Fo] = dynamic_ball_loads(B,ai,ao,Oi,Oo);

dQ = [Qi-Qri;
      Qo-Qro;
      Qi.*sin(ai) - Fi.*cos(ai) - Qo.*sin(ao) + Fo.*cos(ao);
      Qi.*cos(ai) + Fi.*sin(ai) - Qo.*cos(ao) - Fo.*sin(ao) + Fc];