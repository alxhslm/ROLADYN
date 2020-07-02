function [F,V,S] = REB_harris_fast(B,States)

q    = States.qi - States.qo;

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

V = eval_contact_law(B,wons*Oi,wons*Oo,Ar,Az);

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
F.FInt = [];

%stiffnesses
S = struct([]);

function L = eval_contact_law(B,Oi,Oo,Ar0,Az0)
%now find the ball forces
if B.Options.bCentrifugal 
    [Qi,Qo,Xr,Xz,wi,wo] = dynamic_contactlaw_harris(B,Oi,Oo,Ar0,Az0);
    Ar = Ar0 - wo - wi;
    Az = Az0;
    
    [Ai,Ao,ai,ao] = race_geometry(Xz,Xr,Az,Ar);
    dbi = Ai - (B.Geometry.RRacei-B.Geometry.D/2);
    dbo = Ao - (B.Geometry.RRaceo-B.Geometry.D/2);
        
    %dynamic loads
    [Fc,Fi,Fo,Mg] = dynamic_ball_loads(B,ai,ao,Oi,Oo);
else
    if B.Options.bRaceCompliancei || B.Options.bRaceComplianceo
        [wi,wo] = race_compliance_contactlaw_harris(B,Ar0,Az0);

    else
        wi = 0*Ar0;
        wo = 0*Ar0;       
    end
    
    Ar = Ar0 - wi - wo;
    Az = Az0;
        
    A = hypot(Az,Ar);
    
    Xz = (Az./A) .* ((B.Geometry.RRaceo-B.Geometry.D/2) + max(A - B.Geometry.A0,0)/(1 + B.Contact.lambda));
    Xr = (Ar./A) .* ((B.Geometry.RRaceo-B.Geometry.D/2) + max(A - B.Geometry.A0,0)/(1 + B.Contact.lambda));
    
    [Ai,Ao,ai,ao] = race_geometry(Xz,Xr,Az,Ar);
    dbi = Ai - (B.Geometry.RRacei-B.Geometry.D/2);
    dbo = Ao - (B.Geometry.RRaceo-B.Geometry.D/2);
    
         
    dn = hypot(Az,Ar) - B.Geometry.A0;
    Qi = hertz_contactlaw(B.Contact.K,B.Contact.n,dn,B.Contact.tol);
    Qo = Qi;

    Fi = 0*A;
    Fo = 0*A;
    Fc = 0*A;
    Mg = 0*A;
end

L.Qi = Qi;
L.Qo = Qo;

L.Qri = Qi.*cos(ai);
L.Qro = Qo.*cos(ao);

L.dbi = dbi;
L.dbo = dbo;

L.Fi = Fi;
L.Fo = Fo;

L.Az = Az;
L.Ar = Ar;

L.alphai = ai;
L.alphao = ao;

L.Fc = Fc;
L.Mg = Mg;

L.Fri = Qi .* cos(ai)  +  Fi .* sin(ai);
L.Fzi = Qi .* sin(ai)  -  Fi .* cos(ai);
L.Fro = Qo .* cos(ao)  +  Fo .* sin(ao);
L.Fzo = Qo .* sin(ao)  -  Fo .* cos(ao);

L.qi = [permute(L.Fri,[3 4 1 2]);
        permute(L.Fzi,[3 4 1 2])];
    
L.qo = [permute(L.Fro,[3 4 1 2]);
        permute(L.Fzo,[3 4 1 2])];

function L = deriv_contact(B,Oi,Oo,Ar0,Az0)
L = eval_contact_law(B,Oi,Oo,Ar0,Az0);

h = 1E-10;
Lr = eval_contact_law(B,Oi,Oo,Ar0+h,Az0);
Lz = eval_contact_law(B,Oi,Oo,Ar0,Az0+h);

L.Ki = [Lr.qi-L.qi, Lz.qi - L.qi];
L.Ko = [Lr.qo-L.qo, Lz.qo - L.qo];
