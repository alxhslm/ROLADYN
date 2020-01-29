function [F,V,S] = REB_harris_fast(B,States)

q    = States.qi - States.qo;

N = 6;
NPts = size(q,2);

Oi = States.Oi;
Oo = States.Oo;
Ai = States.Ai;
Ao = States.Ao;

E = B.Elements;

%create some vectors of the right size
wons = (E.psi*0+1);
x0 = (q(1,:)*0 + 1);

%if we're not including VC effects, set Acage to 0 to "fix" the cage
if ~B.Options.bVC
    Acage = 0*Ai;
else
    Acage = B.Kinematics.rCagei * Ai + B.Kinematics.rCageo * Ao;
end

wi = 0*wons*x0;
wo = 0*wons*x0;

sgn = sign(E.z/(E.z(1)+eps));
zi = B.Geometry.zRacei*sgn;
Zi = zi*x0;
zo = B.Geometry.zRaceo*sgn;
Zo = zo*x0;

PSI = E.psi*x0 + wons*Acage;
cosPSI = cos(PSI);
sinPSI = sin(PSI);

ALPHA = E.alpha*x0;
cosALPHA = cos(ALPHA);
sinALPHA = sin(ALPHA);

axial = wons*q(3,:);
radial = cos(PSI).*(wons*q(1,:)) + sin(PSI).*(wons*q(2,:));
theta = (sin(PSI).*(wons*q(4,:)) - cos(PSI).*(wons*q(5,:)));

z = axial  + B.Geometry.rRacei * sin(theta) + Zi.*cos(theta) - B.Geometry.cz;
r = radial + B.Geometry.rRacei * cos(theta) - Zi.*sin(theta) - B.Geometry.cr;

dz = z - Zi;
dr = r - B.Geometry.rRacei;

Az = B.Geometry.A0*sinALPHA + dz;
Ar = B.Geometry.A0*cosALPHA + dr;
A = sqrt(Az.^2 + Ar.^2);

Xz0 = (Az./A) .* ((B.Geometry.RRaceo-B.Geometry.D/2) + max(A - B.Geometry.A0,0)/(1 + B.Contact.lambda));
Xr0 = (Ar./A) .* ((B.Geometry.RRaceo-B.Geometry.D/2) + max(A - B.Geometry.A0,0)/(1 + B.Contact.lambda));
[Ai0,Ao0,ai0,ao0] = race_geometry(Xz0,Xr0,Az,Ar);
dn0 = A - B.Geometry.A0;
dbi0 = Ai0 - (B.Geometry.RRacei-B.Geometry.D/2);
dbo0 = Ao0 - (B.Geometry.RRaceo-B.Geometry.D/2);
[Qi0,K0] = hertz_contactlaw(B.Contact.K,B.Contact.n,dn0,B.Contact.tol);
Qo0 = Qi0;

%dynamic loads
[Fc,Fi,Fo,Mg] = dynamic_ball_loads(B,ai0,ao0,wons*Oi,wons*Oo);

%now find the ball forces
if (B.Options.bCentrifugal || B.Options.bGyro) 
    [Qi,Qo,Xr,Xz,wi,wo,Ktoti,Ktoto] = dynamic_contactlaw_harris(B.Contact,B.Geometry,B.Race,B.Options,Fc,Fi,Fo,Ar,Az);
    [Ai,Ao,ai,ao] = race_geometry(Xz,Xr,Az,Ar);
    dbi = Ai - (B.Geometry.RRacei-B.Geometry.D/2) - wi;
    dbo = Ao - (B.Geometry.RRaceo-B.Geometry.D/2) - wo;
    dn = dbi + dbo;
elseif B.Options.bRaceCompliancei || B.Options.bRaceComplianceo
    [Qi,Qo,wi,wo,Ktoti,Ktoto] = race_compliance_contactlaw(B.Contact,B.Race,B.Options,dn0);
    dn = dn0 - (wo + wi);
    dbo = dn/(1+B.Contact.lambda);
    dbi = dn-dbo;
    db = dbo + wo;
    Xz = (B.Geometry.RRaceo-B.Geometry.D/2 + db).*sin(ai0);
    Xr = (B.Geometry.RRaceo-B.Geometry.D/2 + db).*cos(ai0);
    ai = ai0;
    ao = ao0;
    Fi = 0*Fi;
    Fo = 0*Fo;
else
    Xz = Xz0;
    Xr = Xr0;
    dbi = dbi0;
    dbo = dbo0;
    Qi = Qi0;
    Qo = Qo0;
    dn = dn0;
    ai = ai0;
    ao = ao0;
    Ktoti = K0;
    Ktoto = K0;
    Fi = 0*Fi;
    Fo = 0*Fo;
end

%% Introduce scaling factor to account for Sjovall
Qi0 = E.r*Qi0;
Qo0 = E.r*Qo0;

Qi = E.r * Qi;
Qo = E.r * Qo;

Ktoti = E.r * Ktoti;
Ktoto = E.r * Ktoto;

Fc = E.r*Fc;
Fi = E.r*Fi;
Fo = E.r*Fo;
Mg = E.r*Mg;

%% Race loads
cosALPHAI = cos(ai);
sinALPHAI = sin(ai);
cosALPHAO = cos(ao);
sinALPHAO = sin(ao);

cosPSI = permute(cosPSI,[3 4 2 1]);
sinPSI = permute(sinPSI,[3 4 2 1]);

Ki = permute(Ktoti,[3 4 2 1]);
cosALPHAI = permute(cosALPHAI,[3 4 2 1]);
sinALPHAI = permute(sinALPHAI,[3 4 2 1]);
Zi = permute(Zi,[3 4 2 1]);

Ko = permute(Ktoto,[3 4 2 1]);
cosALPHAO = permute(cosALPHAO,[3 4 2 1]);
sinALPHAO = permute(sinALPHAO,[3 4 2 1]);
Zo = permute(Zo,[3 4 2 1]);

Ji = [cosALPHAI.*cosPSI;
      cosALPHAI.*sinPSI;
      sinALPHAI;
    - Zi.*cosALPHAI.*sinPSI + B.Geometry.rRacei*sinALPHAI.*sinPSI
      Zi.*cosALPHAI.*cosPSI - B.Geometry.rRacei*sinALPHAI.*cosPSI;
      0*Zi];

Jo =-[cosALPHAO.*cosPSI;
      cosALPHAO.*sinPSI;
      sinALPHAO;
    - Zo.*cosALPHAO.*sinPSI + B.Geometry.rRaceo*sinALPHAO.*sinPSI
      Zo.*cosALPHAO.*cosPSI - B.Geometry.rRaceo*sinALPHAO.*cosPSI;
      0*Zo];

Wi =  permute(sum(Ji.*permute(Qi,[3 4 2 1]),4),[1 3 2]);
Wo =  permute(sum(Jo.*permute(Qo,[3 4 2 1]),4),[1 3 2]);

%forces
F.Fi = Wi;
F.Fo = Wo;
F.FInt = [];

%contact loads
V.Qi = Qi; V.Qo = Qo; 
V.Fi = Fi; V.Fo = Fo; 
V.Qi0 = Qi0; V.Qo0 = Qo0;

%geometry
V.alpha_i = ai; V.alpha_o = ao;
V.dbi = dbi; V.dbo = dbo;
V.dn = dn;
V.Xr = Xr; V.Xz = Xz;
V.Ar = Ar; V.Az = Az;

%dynamic loads
V.Fc = Fc;   
V.Mg = Mg;

%race compliance
V.wi  = wi;  V.wo = wo;

%contact stiffness
V.Ki = Ktoti; V.Ko = Ktoto;

%stiffnesses
if nargout > 2
    S.Kqiqi = sum(Ji.*Ki.*permute(Ji,[2 1 3 4]),4); 
    S.Kqoqi = sum(Jo.*Ko.*permute(Ji,[2 1 3 4]),4);
    S.Kqiqo = sum(Ji.*Ki.*permute(Jo,[2 1 3 4]),4);
    S.Kqoqo = sum(Jo.*Ko.*permute(Jo,[2 1 3 4]),4);
    
    S.Kqiqi(:,[4 5],:) = 0;
    S.Kqoqi(:,[4 5],:) = 0;
    S.Kqiqo(:,[4 5],:) = 0;
    S.Kqoqo(:,[4 5],:) = 0;
    
    S.Kqiqi([4 5],:,:) = 0;
    S.Kqoqi([4 5],:,:) = 0;
    S.Kqiqo([4 5],:,:) = 0;
    S.Kqoqo([4 5],:,:) = 0;

    S.Kxqi = zeros(0,N,NPts);
    S.Kxqo = zeros(0,N,NPts);
    S.Kqix = zeros(N,0,NPts);
    S.Kqox = zeros(N,0,NPts);
    S.Kxx  = [];

    %damping
    S.Cqiqi = zeros(N,N,NPts);
    S.Cqoqi = zeros(N,N,NPts);
    S.Cqiqo = zeros(N,N,NPts);
    S.Cqoqo = zeros(N,N,NPts);

    S.Cxqi = zeros(0,N,NPts);
    S.Cxqo = zeros(0,N,NPts);
    S.Cqix = zeros(N,0,NPts);
    S.Cqox = zeros(N,0,NPts);
    S.Cxx  = [];
end