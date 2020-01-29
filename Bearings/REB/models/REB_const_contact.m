function [F,V,S] = REB_const_contact(B,States)

q    = States.qi - States.qo;
xInt = States.xInt;

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

vr = 0*wons*x0;
wi = 0*wons*x0;
wo = 0*wons*x0;

count = 0;
if B.Options.bCentrifugal
    vr = xInt(count + (1:E.N),:);
    count  = count + E.N;
end
if B.Options.bRaceCompliancei
    wi = xInt(count+(1:E.N),:);
    count  = count + E.N;
end
if B.Options.bRaceComplianceo
    wo = xInt(count+(1:E.N),:);
    count  = count + E.N;
end

sgn = sign(E.z/(E.z(1)+eps));
z = B.Geometry.zRacei*sgn;
Z = z*x0;

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

dz = z - Z;
dr = r - B.Geometry.rRacei;

dn0 = dr .* cosALPHA + dz .* sinALPHA;
db0 = dn0  / (1 + B.Contact.lambda);
dbi0 = dn0 - db0;
dbo0 = db0;
Qi0 = hertz_contactlaw(B.Contact.K,B.Contact.n,dn0,B.Contact.tol);
Qo0 = Qi0;

%contact angles are constant by definition
ai = E.alpha*x0;
ao = E.alpha*x0;

%now find the ball forces
if B.Options.bCentrifugal
    db = vr;
    dbi = dn0-db-wi;
    dbo = db-wo;
    Qi = hertz_contactlaw(B.Contact.Inner.K,B.Contact.n,dbi,B.Contact.tol);
    Qo = hertz_contactlaw(B.Contact.Outer.K,B.Contact.n,dbo,B.Contact.tol);
    dn = dbi + dbo;
elseif B.Options.bRaceCompliancei || B.Options.bRaceComplianceo
    dn =  (dn0 - (wo + wi));
    Qi = hertz_contactlaw(B.Contact.K,B.Contact.n,dn,B.Contact.tol);
    Qo = Qi;
    db = dn/(1+B.Contact.lambda);
    dbo = db;
    dbi = dn - db;   
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

%dynamic loads
Fc = dynamic_ball_loads(B,ai,ao,wons*Oi,wons*Oo);
if B.Options.bCentrifugal 
    fErr = Qi - Qo + Fc./cos(ao);
else
    fErr = [];
end

% race compliance
if B.Options.bRaceCompliancei
    Qri = -race_compliance_loads(B.Race.Inner,-wi);
    fErr = [fErr;
            Qi - Qri];
end
if B.Options.bRaceComplianceo
    Qro = race_compliance_loads(B.Race.Outer, wo);
    fErr = [fErr;
            Qo - Qro];
end

%% Introduce scaling factor to account for Sjovall
Qi0 = E.r*Qi0;
Qo0 = E.r*Qo0;

Qi = E.r * Qi;
Qo = E.r * Qo;

Fc = E.r*Fc;

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
F.FInt = fErr;

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