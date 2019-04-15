function [F,V,S] = REB_const_contact_fast(B,States)

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
z = B.Geometry.zRacei*sgn;
Z = z*x0;

PSI = E.psi*x0 + wons*Acage;
cosPSI = cos(PSI);
sinPSI = sin(PSI);

ALPHA = E.alpha*x0;
cosALPHA = cos(ALPHA);
sinALPHA = sin(ALPHA);

dz  = wons*q(3,:) + B.Geometry.rRacei*(sinPSI.*(wons*q(4,:)) - cosPSI.*(wons*q(5,:))) - B.Geometry.cz;
dr  = cosPSI.*(wons*q(1,:)) + sinPSI.*(wons*q(2,:)) - Z.*(sinPSI.*(wons*q(4,:)) - cosPSI.*(wons*q(5,:))) - B.Geometry.cr;

dn0 = dr .* cosALPHA + dz .* sinALPHA;
db0 = dn0  / (1 + B.Contact.lambda);
dbi0 = dn0 - db0;
dbo0 = db0;

[Qi0,K0] = hertz_contactlaw(B.Contact.K,B.Contact.n,dn0,B.Contact.tol);
Qo0 = Qi0;

%contact angles are constant by definition
ai = E.alpha*x0;
ao = E.alpha*x0;

%dynamic loads
Fc = dynamic_ball_loads(B,ai,ao,wons*Oi,wons*Oo);

%now find the ball forces
if B.Options.bCentrifugal
    [Qi,Qo,vr,wi,wo,Ktot] = dynamic_contactlaw(B.Contact,B.Race,B.Options,(Fc./cosALPHA),dn0);
    db = vr;
    dbi = dn0-db-wi;
    dbo = db-wo;
    dn = dn0 - (wo + wi);
elseif B.Options.bRaceCompliancei || B.Options.bRaceComplianceo
    [Qi,Qo,wi,wo,Ktot] = race_compliance_contactlaw(B.Contact,B.Race,B.Options,dn0);
    dn = dn0 - (wo + wi);
    dbo = dn/(1+B.Contact.lambda);
    dbi = dn - dbo;
    db = dbo + wo;
else
    db = db0;
    dn = dn0;
    dbi = dbi0;
    dbo = dbo0;
    Qi = Qi0;
    Qo = Qo0;
    Ktot = K0;
end

%% Introduce scaling factor to account for Sjovall
Qi0 = E.r*Qi0;
Qo0 = E.r*Qo0;

Qi = E.r * Qi;
Qo = E.r * Qo;

Ktot = E.r * Ktot;

Fc = E.r*Fc;

%% Geometric parameters

Az = B.Geometry.A0*sinALPHA + dn.*sinALPHA;
Ar = B.Geometry.A0*cosALPHA + dn.*cosALPHA;
Xz = (B.Geometry.RRaceo-B.Geometry.D/2)*sinALPHA + db.*sinALPHA;
Xr = (B.Geometry.RRaceo-B.Geometry.D/2)*cosALPHA + db.*cosALPHA;

%% Race loads
Ki = permute(Ktot,[3 4 2 1]);
cosPSI = permute(cosPSI,[3 4 2 1]);
sinPSI = permute(sinPSI,[3 4 2 1]);
cosALPHA = permute(cosALPHA,[3 4 2 1]);
sinALPHA = permute(sinALPHA,[3 4 2 1]);
Z = permute(Z,[3 4 2 1]);

J = [cosALPHA.*cosPSI;
     cosALPHA.*sinPSI;
     sinALPHA;
    - Z.*cosALPHA.*sinPSI + B.Geometry.rRacei*sinALPHA.*sinPSI
      Z.*cosALPHA.*cosPSI - B.Geometry.rRacei*sinALPHA.*cosPSI;
      0*Z];

Wi =  permute(sum(J.*permute(Qi,[3 4 2 1]),4),[1 3 2]);
Wo = -permute(sum(J.*permute(Qo,[3 4 2 1]),4),[1 3 2]);

%forces
F.Fi = Wi;
F.Fo = Wo;
F.FInt = [];

%contact loads
V.Qi = Qi; V.Qo = Qo; 
V.Qi0 = Qi0; V.Qo0 = Qo0;

%geometry
V.alpha_i = ai; V.alpha_o = ao;
V.dbi = dbi; V.dbo = dbo;
V.dn = dn;
V.Xr = Xr; V.Xz = Xz;
V.Ar = Ar; V.Az = Az;

%dynamic loads
V.Fc = Fc;   

%race compliance
V.wi  = wi;  V.wo = wo;

%contact stiffness
V.Ki = Ktot; V.Ko = Ktot;

%stiffnesses
if nargout > 2
    K = sum(J.*Ki.*permute(J,[2 1 3 4]),4);

    S.Kqiqi = K;
    S.Kqoqi = K;
    S.Kqiqo = K;
    S.Kqoqo = K;

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