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
Ar = B.Geometry.A0*cosALPHA + dr - wi - wo;
A = sqrt(Az.^2 + Ar.^2);

%now find the ball forces
if B.Options.bCentrifugal
    Xz = (B.Geometry.RRaceo-B.Geometry.D/2)*sinALPHA*x0 - vz.*(sgn*x0);
    Xr = (B.Geometry.RRaceo-B.Geometry.D/2)*cosALPHA*x0 + vr;
    
    [Ai,Ao,V.alphai,V.alphao] = race_geometry(Xz,Xr,Az,Ar);
       
    V.dbi = Ai - (B.Geometry.RRacei-B.Geometry.D/2);
    V.dbo = Ao - (B.Geometry.RRaceo-B.Geometry.D/2);
    V.Qi = hertz_contactlaw(B.Contact.Inner.K,B.Contact.n,V.dbi,B.Contact.tol);
    V.Qo = hertz_contactlaw(B.Contact.Outer.K,B.Contact.n,V.dbo,B.Contact.tol);
else
    Xz = (Az./A) .* ((B.Geometry.RRaceo-B.Geometry.D/2) + max(A - B.Geometry.A0,0)/(1 + B.Contact.lambda));
    Xr = (Ar./A) .* ((B.Geometry.RRaceo-B.Geometry.D/2) + max(A - B.Geometry.A0,0)/(1 + B.Contact.lambda));
    [V.Ai,V.Ao,V.alphai,V.alphao] = race_geometry(Xz,Xr,Az,Ar);
    dn = A-B.Geometry.A0;
    V.dbi = Ai - (B.Geometry.RRacei-B.Geometry.D/2);
    V.dbo = Ao - (B.Geometry.RRaceo-B.Geometry.D/2);
    V.Qi = hertz_contactlaw(B.Contact.K,B.Contact.n,dn,B.Contact.tol);
    V.Qo = V.Qi;
end

%dynamic loads
[V.Fc,V.Fi,V.Fo,V.Mg] = dynamic_ball_loads(B,V.alphai,V.alphao,wons*Oi,wons*Oo);

if B.Options.bCentrifugal
    fErr = [V.Qi.*sin(V.alphai) - V.Fi.*cos(V.alphai) - V.Qo.*sin(V.alphao) + V.Fo.*cos(V.alphao);
            V.Qi.*cos(V.alphai) + V.Fi.*sin(V.alphai) - V.Qo.*cos(V.alphao) - V.Fo.*sin(V.alphao) + V.Fc];
else
    fErr = [];
end

%race compliance
if B.Options.bRaceCompliancei
    V.Qri = race_compliance_loads(B.Race.Inner,-wi);
    fErr = [fErr;
            V.Qi.*cos(V.alphai) - V.Qri];
else
    V.Qri = V.Qi.*cos(V.alphai);
end
if B.Options.bRaceComplianceo
    V.Qro = race_compliance_loads(B.Race.Outer, wo);
    fErr = [fErr;
            V.Qo.*cos(V.alphao) - V.Qro];
else
    V.Qro = V.Qo.*cos(V.alphao);
end

%% Introduce scaling factor to account for Sjovall
V.Qi = E.r * V.Qi;
V.Qo = E.r * V.Qo;

V.Qri = E.r * V.Qri;
V.Qro = E.r * V.Qro;

V.Fc = E.r*V.Fc;
V.Fi = E.r*V.Fi;
V.Fo = E.r*V.Fo;
V.Mg = E.r*V.Mg;

%% Race loads
Fri = V.Qi.*cos(V.alphai) + V.Fi.*sin(V.alphai);
Fzi = V.Qi.*sin(V.alphai) - V.Fi.*cos(V.alphai);
Wi = [sum(Fri.*cosPSI);
      sum(Fri.*sinPSI);
      sum(Fzi);
      sum( B.Geometry.rRacei.*Fzi.*sinPSI - Zi.*Fri.*sinPSI);
      sum(-B.Geometry.rRacei.*Fzi.*cosPSI + Zi.*Fri.*cosPSI);
      0*x0];

Fro = V.Qo.*cos(V.alphao) + V.Fo.*sin(V.alphao);
Fzo = V.Qo.*sin(V.alphao) - V.Fo.*cos(V.alphao);
Wo =-[sum(Fro.*cosPSI);
      sum(Fro.*sinPSI);
      sum(Fzo);
      sum( B.Geometry.rRaceo.*Fzo.*sinPSI - Zo.*Fro.*sinPSI);
      sum(-B.Geometry.rRaceo.*Fzo.*cosPSI + Zo.*Fro.*cosPSI);
      0*x0];
  
%forces
F.Fi = Wi;
F.Fo = Wo;
F.FInt = fErr;
 
%stiffnesses
S = struct([]);