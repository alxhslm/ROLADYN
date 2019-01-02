function [F,V,S] = REB_const_contact_fast(B,States)

q    = States.qi    - States.qo;
xInt = States.xInt;

% if B.bPinned
%     q(4:5,:) = 0;
% end

Oi = States.Oi;
Oo = States.Oo;
Ai = States.Ai;
Ao = States.Ao;

E = B.Elements;

%create some vectors of the right size
wons = (E.psi*0+1);
x0 = (q(1,:)*0 + 1);

%if we're not including VC effects, set Ocage to 0 to "fix" the cage
if ~B.Options.bVC
    Acage = 0*Ai;
else
    Acage = B.Kinematics.rCagei * Ai + B.Kinematics.rCageo * Ao;
end
Ocage = B.Kinematics.rCagei * Oi + B.Kinematics.rCageo * Oo;

vr = 0*wons*x0;
wi = 0*wons*x0;
wo = 0*wons*x0;

sgn = sign(E.z/(E.z(1)+eps));
z = B.Geometry.zi*sgn;
Z = z*x0;

PSI = E.psi*x0 + Acage;

dz  = wons*q(3,:) + B.Geometry.Ri*(sin(PSI).*(wons*q(4,:)) - cos(PSI).*(wons*q(5,:))) - B.Geometry.cz;
dr  = cos(PSI).*(wons*q(1,:)) + sin(PSI).*(wons*q(2,:)) - Z.*(sin(PSI).*(wons*q(4,:)) - cos(PSI).*(wons*q(5,:))) - B.Geometry.cr;

dn = dr .* (cos(E.alpha)*x0) + dz .* (sin(E.alpha)*x0);
lambda = (B.Contact.Outer.K / B.Contact.Inner.K)^(1/B.Contact.n);
db0 = dn / (1 + lambda);
dbi0 = dn-db0;
dbo0 = db0;

tol = 1E-8;
Qi0 = E.r*hertz_contact(B.Contact.K,B.Contact.n,dn,tol);
Qo0 = Qi0;

%contact angles are constant by definition
alpha_i = E.alpha*x0;
alpha_o = E.alpha*x0;

%dynamic loads
Fc = E.r*dynamic_ball_loads(B,alpha_i,alpha_o,wons*Oi,wons*Oo);

%now find the ball forces
if B.Options.bCentrifugal
    db = vr;
    dbi = dn-db-wi;
    dbo = db-wo;
    Qi = E.r * dynamic_contact(B.Contact,Fc/E.r,dn,tol);
    Qo = Qi + Fc;
elseif B.Options.bRaceCompliancei || B.Options.bRaceComplianceo
    Qi = E.r*race_contact(B.Contact,B.Race,B.Options,dn,tol);
    Qo = Qi;
    db = dn/(1+lambda);
    dbo = db;
    dbi = dn - db;   
else
    db = db0;
    
    dbi = dbi0;
    dbo = dbo0;
    
    Qi = Qi0;
    Qo = Qo0;
end

Fi = 0*Qi;
Fo = 0*Qo;

%race compliance
Qri = Qi;
Qro = Qo;

Xz = B.Geometry.D/2*sin(E.alpha)*x0 + db.*sin(alpha_i);
Xr = B.Geometry.D/2*cos(E.alpha)*x0 + db.*cos(alpha_i);

%% Race loads
Fri = Qi.*cos(alpha_i) + Fi.*sin(alpha_i);
Fzi = Qi.*sin(alpha_i) - Fi.*cos(alpha_i);
Wi = [sum(Fri.*cos(PSI));
    sum(Fri.*sin(PSI));
    sum(Fzi);
    sum( B.Geometry.Ri.*Fzi.*sin(PSI) - Z.*Fri.*sin(PSI));
    sum(-B.Geometry.Ri.*Fzi.*cos(PSI) + Z.*Fri.*cos(PSI));
    0*x0];

Fro = Qo.*cos(alpha_o) + Fo.*sin(alpha_o);
Fzo = Qo.*sin(alpha_o) - Fo.*cos(alpha_o);
Wo =-[sum(Fro.*cos(PSI));
      sum(Fro.*sin(PSI));
      sum(Fzo);
      sum( B.Geometry.Ro.*Fzo.*sin(PSI) - Z.*Fro.*sin(PSI));
      sum(-B.Geometry.Ro.*Fzo.*cos(PSI) + Z.*Fro.*cos(PSI));
      0*x0];

%forces
F.Fi = Wi;
F.Fo = Wo;
F.FInt = [];

%contact loads
V.Qi = Qi; V.Qo = Qo; 
V.Fi = Fi; V.Fo = Fo;

%geometry
V.alpha_i = alpha_i; V.alpha_o = alpha_o;
V.dbi = dbi; V.dbo = dbo;
V.dn = dn;
V.Ai = Ai;  V.Ao = Ao;
V.Xr = Xr;  
V.Xz = Xz;

%dynamic loads
V.Fc = Fc;   

%race compliance
V.Qri = Qri; V.Qro = Qro;
V.wi  = wi;  V.wo = wo;
 
%stiffnesses
S = struct([]);