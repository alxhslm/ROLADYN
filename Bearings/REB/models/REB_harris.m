function [F,V,S] = REB_harris(B,States)

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

vz = 0*wons*x0;
vr = 0*wons*x0;
wi = 0*wons*x0;
wo = 0*wons*x0;

count = 0;
if B.Options.bCentrifugal
    vz = xInt(count + (1:E.N),:);
    vr = xInt(count + E.N + (1:E.N),:);
    count  = count + 2*E.N;
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
Z = B.Geometry.zRacei*sgn*x0;

PSI = E.psi*x0 + wons*Acage;
cosPSI = cos(PSI);
sinPSI = sin(PSI);

ALPHA = E.alpha*x0;
cosALPHA = cos(ALPHA);
sinALPHA = sin(ALPHA);

axial = wons*q(3,:);
radial = cos(PSI).*(wons*q(1,:)) + sin(PSI).*(wons*q(2,:));
theta = (sin(PSI).*(wons*q(4,:)) - cos(PSI).*(wons*q(5,:)));

z = axial + B.Geometry.rRacei * sin(theta) + Z.*cos(theta)  - B.Geometry.cz;
r = radial + B.Geometry.rRacei * cos(theta) - Z.*sin(theta) - B.Geometry.cr;
% dz  = axial + B.Geometry.rRacei*theta - B.Geometry.cz;
% dr  = radial - Z.*theta - B.Geometry.cr;

dz = z - Z;
dr = r - B.Geometry.rRacei;

Az = B.Geometry.A0*sinALPHA + dz;
Ar = B.Geometry.A0*cosALPHA + dr;
A = sqrt(Az.^2 + Ar.^2);

Xz0 = (Az./A) .* ((B.Geometry.RRaceo-B.Geometry.D/2) + max(A - B.Geometry.A0,0)/(1 + B.Contact.lambda));
Xr0 = (Ar./A) .* ((B.Geometry.RRaceo-B.Geometry.D/2) + max(A - B.Geometry.A0,0)/(1 + B.Contact.lambda));
[Ai0,Ao0,alpha_i0,alpha_o0] = race_geometry(Xz0,Xr0,Az,Ar);
dbi0 = Ai0 - (B.Geometry.RRacei-B.Geometry.D/2) - wi;
dbo0 = Ao0 - (B.Geometry.RRaceo-B.Geometry.D/2) - wo;
db0 = A-B.Geometry.A0-wi-wo;
Qi0 = hertz_contact(E.r*B.Contact.K,B.Contact.n,db0,B.Contact.tol);
Qo0 = Qi0;

%now find the ball forces
if (B.Options.bCentrifugal || B.Options.bGyro) 
    Xz = (B.Geometry.RRaceo-B.Geometry.D/2)*sinALPHA*x0 - vz.*(sgn*x0);
    Xr = (B.Geometry.RRaceo-B.Geometry.D/2)*cosALPHA*x0 + vr;
    
    [Ai,Ao,alpha_i,alpha_o] = race_geometry(Xz,Xr,Az,Ar);
       
    dbi = Ai - (B.Geometry.RRacei-B.Geometry.D/2) - wi;
    dbo = Ao - (B.Geometry.RRaceo-B.Geometry.D/2) - wo;
    Qi = hertz_contactlaw(E.r*B.Contact.Inner.K,B.Contact.n,dbi,B.Contact.tol);
    Qo = hertz_contactlaw(E.r*B.Contact.Outer.K,B.Contact.n,dbo,B.Contact.tol);

    dn = dbi + dbo;
else
    alpha_i = alpha_i0;
    alpha_o = alpha_o0;

    if B.Options.bRaceCompliancei || B.Options.bRaceComplianceo
        dn =  (dn0 - (wo + wi));
        Qi = E.r*hertz_contactlaw(B.Contact.K,B.Contact.n,dn,B.Contact.tol);
        Qo = Qi;
        db = dn/(1+lambda);
        dbo = db;
        dbi = dn - db;  

        Xz = (B.Geometry.RRaceo-B.Geometry.D/2 + db)*sin(alpha_i0);
        Xr = (B.Geometry.RRaceo-B.Geometry.D/2 + db)*cos(alpha_i0);

        dn = dbi + dbo;
    else
        Xz = Xz0;
        Xr = Xr0;
        dbi = dbi0;
        dbo = dbo0;
        Qi = Qi0;
        Qo = Qo0;
        dn = dn0;
    end
end

%dynamic loads
[Fc,Fi,Fo,Mg] = dynamic_ball_loads(B,alpha_i,alpha_o,wons*Oi,wons*Oo);
Fc = E.r*Fc; Fi = E.r*Fi; Fo = E.r*Fo; Mg = E.r*Mg;

if (B.Options.bCentrifugal || B.Options.bGyro)
    fErr = [Qi.*sin(alpha_i) - Fi.*cos(alpha_i) - Qo.*sin(alpha_o) + Fo.*cos(alpha_o);
            Qi.*cos(alpha_i) + Fi.*sin(alpha_i) - Qo.*cos(alpha_o) - Fo.*sin(alpha_o) + Fc];
else
    fErr = [];
end

if ~B.Options.bGyro
    Fi = 0*Fi;
    Fo = 0*Fo;
end

%race compliance
if B.Options.bRaceCompliancei
    Qri = race_compliance_loads(B.Race.Inner,-wi);
    fErr = [fErr;
            Qi.*cos(alpha_i) - Qri];
end
if B.Options.bRaceComplianceo
    Qro = race_compliance_loads(B.Race.Outer, wo);
    fErr = [fErr;
            Qo.*cos(alpha_o) - Qro];
end

%% Race loads
Fri = Qi.*cos(alpha_i) + Fi.*sin(alpha_i);
Fzi = Qi.*sin(alpha_i) - Fi.*cos(alpha_i);
Wi = [sum(Fri.*cosPSI);
      sum(Fri.*sinPSI);
      sum(Fzi);
      sum( B.Geometry.rRacei.*Fzi.*sinPSI - Z.*Fri.*sinPSI);
      sum(-B.Geometry.rRacei.*Fzi.*cosPSI + Z.*Fri.*cosPSI);
      0*x0];

Fro = Qo.*cos(alpha_o) + Fo.*sin(alpha_o);
Fzo = Qo.*sin(alpha_o) - Fo.*cos(alpha_o);
Wo =-[sum(Fro.*cosPSI);
      sum(Fro.*sinPSI);
      sum(Fzo);
      sum( B.Geometry.rRaceo.*Fzo.*sinPSI - Z.*Fro.*sinPSI);
      sum(-B.Geometry.rRaceo.*Fzo.*cosPSI + Z.*Fro.*cosPSI);
      0*x0];

%forces
F.Fi = Wi;
F.Fo = Wo;
F.FInt = fErr;

%contact loads
V.Qi = Qi; V.Qo = Qo; 
V.Fi = Fi; V.Fo = Fo;
V.Qi0 = Qi0; V.Qo0 = Qo0;

%geometry
V.alpha_i = alpha_i; V.alpha_o = alpha_o;
V.dbi = dbi; V.dbo = dbo;
V.dn = dn;
V.Ar = Ar; V.Az = Az;
V.Xr = Xr; V.Xz = Xz;

V.alpha_i0 = alpha_i0; V.alpha_o0 = alpha_o0;
V.dbi0 = dbi0; V.dbo0 = dbo0;
V.Xr0 = Xr0;  V.Xz0 = Xz0;

%dynamic loads
V.Fc = Fc; V.Mg = Mg;

%race compliance
V.wi  = wi;  V.wo = wo;
 
%stiffnesses
S = struct([]);