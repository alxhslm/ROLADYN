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
Qi0 = E.r*hertz_contact(B.Contact.K,B.Contact.n,dn0,B.Contact.tol);
Qo0 = Qi0;

%contact angles are constant by definition
ai = E.alpha*x0;
ao = E.alpha*x0;

%dynamic loads
Fc = E.r*dynamic_ball_loads(B,ai,ao,wons*Oi,wons*Oo);

%now find the ball forces
%TODO: this currently includes EITHER centrifugal loads OR race compliance. There needs to be an option to include both.
if B.Options.bCentrifugal
    [Qi,Qo,vr] = dynamic_contact(B.Contact,Fc/E.r,dn0,tol);
    Qi = E.r * Qi;
    Qo = E.r * Qo;
    db = vr;
    dbi = dn0-db;
    dbo = db;
    dn = dbi + dbo;
elseif B.Options.bRaceCompliancei || B.Options.bRaceComplianceo
    [Qi,Qo,wi,wo] = race_contact(B.Contact,B.Race,B.Options,dn0,tol);
    Qi = E.r*Qi;
    Qo = E.r*Qo;
    dbo = (dn0-wi-wo)/(1+lambda);
    db = dbo + wo;
    dbi = dn0-db-wi;
    dn = dn0 - (wo + wi);
else
    db = db0;
    dn = dn0;
    dbi = dbi0;
    dbo = dbo0;
    Qi = Qi0;
    Qo = Qo0;
end

Az = B.Geometry.A0*sinALPHA + dn.*sinALPHA;
Ar = B.Geometry.A0*cosALPHA + dn.*cosALPHA;
Xz = (B.Geometry.RRaceo-B.Geometry.D/2)*sinALPHA + db.*sinALPHA;
Xr = (B.Geometry.RRaceo-B.Geometry.D/2)*cosALPHA + db.*cosALPHA;

%% Race loads
Wi = [sum(Qi.*cosALPHA.*cosPSI);
      sum(Qi.*cosALPHA.*sinPSI);
      sum( Qi.*sinALPHA);
      sum( B.Geometry.rRacei.*Qi.*sinALPHA.*sinPSI - Z.*Qi.*cosALPHA.*sinPSI);
      sum(-B.Geometry.rRacei.*Qi.*sinALPHA.*cosPSI + Z.*Qi.*cosALPHA.*cosPSI);
      0*x0];

Wo =-[sum(Qo.*cosALPHA.*cosPSI);
      sum(Qo.*cosALPHA.*sinPSI);
      sum(Qo.*sinALPHA);
      sum( B.Geometry.rRaceo.*Qo.*sinALPHA.*sinPSI - Z.*Qo.*cosALPHA.*sinPSI);
      sum(-B.Geometry.rRaceo.*Qo.*sinALPHA.*cosPSI + Z.*Qo.*cosALPHA.*cosPSI);
      0*x0];

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
S = struct([]);