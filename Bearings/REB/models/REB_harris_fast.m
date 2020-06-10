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

V = deriv_contact(B,Oi*wons,Oo*wons,Ar,Az);

%% Introduce scaling factor to account for Sjovall
V.Qi = E.r * V.Qi;
V.Qo = E.r * V.Qo;

V.Ki = E.r * V.Ki;
V.Ko = E.r * V.Ko;

V.Fc = E.r*V.Fc;
V.Fi = E.r*V.Fi;
V.Fo = E.r*V.Fo;
V.Mg = E.r*V.Mg;

%% Race loads
cosPSI = permute(cosPSI,[3 4 2 1]);
sinPSI = permute(sinPSI,[3 4 2 1]);

Zi = permute(Zi,[3 4 2 1]);
Zo = permute(Zo,[3 4 2 1]);

Ji = [cosPSI 0*Zi;
      sinPSI 0*Zi;
      0*Zi   0*Zi+1;
    - Zi.*sinPSI, + B.Geometry.rRacei.*sinPSI
      Zi.*cosPSI, - B.Geometry.rRacei.*cosPSI;
      0*Zi 0*Zi];

Jo =-[cosPSI 0*Zo;
      sinPSI 0*Zo;
      0*Zo   0*Zo+1;
    - Zo.*sinPSI,   B.Geometry.rRaceo.*sinPSI
      Zo.*cosPSI, - B.Geometry.rRaceo.*cosPSI;
      0*Zo  0*Zi];

Wi =  permute(sum(mtimesx(Ji,permute(V.qi,[1 2 4 3])),4),[1 3 2]);
Wo =  permute(sum(mtimesx(Jo,permute(V.qo,[1 2 4 3])),4),[1 3 2]);

%forces
F.Fi = Wi;
F.Fo = Wo;
F.FInt = [];

%stiffnesses
if nargout > 2
    Jit = permute(Ji,[2 1 3 4]);
    Jot = permute(Jo,[2 1 3 4]);
    
    S.Kqiqi = sum(mtimesx(Ji,mtimesx(permute(V.Ki,[1 2 4 3]),Jit)),4); 
    S.Kqoqi = sum(mtimesx(Jo,mtimesx(permute(V.Ko,[1 2 4 3]),Jit)),4); 
    S.Kqiqo = sum(mtimesx(Ji,mtimesx(permute(V.Ki,[1 2 4 3]),Jot)),4); 
    S.Kqoqo = sum(mtimesx(Jo,mtimesx(permute(V.Ko,[1 2 4 3]),Jot)),4); 
    
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


function L = eval_contact_law(B,Oi,Oo,Ar0,Az0)
%now find the ball forces
if (B.Options.bCentrifugal || B.Options.bGyro) 
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

L.dbi = dbi;
L.dbo = dbo;

L.Fi = Fi;
L.Fo = Fo;

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