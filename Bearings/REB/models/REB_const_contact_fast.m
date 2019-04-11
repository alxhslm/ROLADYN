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
lambda = (B.Contact.Outer.K / B.Contact.Inner.K)^(1/B.Contact.n);
db0 = dn0  / (1 + lambda);
dbi0 = dn0 - db0;
dbo0 = db0;

[Qi0,K0] = hertz_contact(B.Contact.K,B.Contact.n,dn0,B.Contact.tol);
Qi0 = E.r*Qi0;
K0 = E.r*K0;
Qo0 = Qi0;

%contact angles are constant by definition
ai = E.alpha*x0;
ao = E.alpha*x0;

%dynamic loads
Fc = E.r*dynamic_ball_loads(B,ai,ao,wons*Oi,wons*Oo);

%now find the ball forces
%TODO: this currently includes EITHER centrifugal loads OR race compliance. There needs to be an option to include both.
if B.Options.bCentrifugal
    [Qi,Qo,vr,Ki,Ko] = dynamic_contact(B.Contact,Fc/E.r,dn0,tol);
    Qi = E.r * Qi;
    Qo = E.r * Qo;
    Ki = E.r * Ki;
    Ko = E.r * Ko;
    
    db = vr;
    dbi = dn0-db;
    dbo = db;
    dn = dbi + dbo;

    Ktot = 1./(1./Ki + 1./Ko);
elseif B.Options.bRaceCompliancei || B.Options.bRaceComplianceo
    [Qi,Qo,wi,wo,Kri,Kro] = race_contact(B.Contact,B.Race,B.Options,dn0,tol);
    Qi = E.r*Qi;
    Qo = E.r*Qo;
    dbo = (dn0-wi-wo)/(1+lambda);
    db = dbo + wo;
    dbi = dn0-db-wi;
    dn = dn0 - (wo + wi);
    Ktot = 1./(1./K0 + 1./Kri + 1./Kro);
else
    db = db0;
    dn = dn0;
    dbi = dbi0;
    dbo = dbo0;
    Qi = Qi0;
    Qo = Qo0;
    Ktot = K0;
end

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