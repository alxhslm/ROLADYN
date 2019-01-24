function [F,V,S] = REB_const_contact(B,States)

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

PSI = E.psi*x0 + Acage;

dz  = wons*q(3,:) + B.Geometry.rRacei*(sin(PSI).*(wons*q(4,:)) - cos(PSI).*(wons*q(5,:))) - B.Geometry.cz;
dr  = cos(PSI).*(wons*q(1,:)) + sin(PSI).*(wons*q(2,:)) - Z.*(sin(PSI).*(wons*q(4,:)) - cos(PSI).*(wons*q(5,:))) - B.Geometry.cr;

dn = dr .* (cos(E.alpha)*x0) + dz .* (sin(E.alpha)*x0);

if B.Options.bSaturate
    rContact = dn>0;
    
    if B.Options.bCentrifugal
        Fc = (B.Dynamics.mb * Ocage.^2 * B.Geometry.dm/2);
        xc =  wons*(Fc/B.Contact.Outer.K).^(1/B.Contact.Outer.n);
        vr = vr.*rContact + (1-rContact).*xc;
    end
    
    wi = wi.*rContact;
    wo = wo.*rContact;
end

lambda = (B.Contact.Outer.K / B.Contact.Inner.K)^(1/B.Contact.n);
db0 = dn / (1 + lambda);
dbi0 = dn-db0;
dbo0 = db0;
Qi0 = E.r*hertz_contact(B.Contact.K,B.Contact.n,dn,0);
Qo0 = Qi0;

%contact angles are constant by definition
alpha_i = E.alpha*x0;
alpha_o = E.alpha*x0;

%now find the ball forces
tol = 1E-8;
if B.Options.bCentrifugal
    db = vr;
    dbi = dn-db-wi;
    dbo = db-wo;
    Qi = E.r*hertz_contact(B.Contact.Inner.K,B.Contact.n,dbi,tol);
    Qo = E.r*hertz_contact(B.Contact.Outer.K,B.Contact.n,dbo,tol);
elseif B.Options.bRaceCompliancei || B.Options.bRaceComplianceo
    dn_race =  (dn - (wo + wi));
    Qi = E.r*hertz_contact(B.Contact.K,B.Contact.n,dn_race,tol);
    Qo = Qi;
    db = dn_race/(1+lambda);
    dbo = db;
    dbi = dn_race - db;   
else
    db = db0;
    
    dbi = dbi0;
    dbo = dbo0;
    
    Qi = Qi0;
    Qo = Qo0;
end

%dynamic loads
if B.Options.bCentrifugal 
    Fc = E.r*dynamic_ball_loads(B,alpha_i,alpha_o,wons*Oi,wons*Oo);
    fErr = Qi - Qo + Fc./cos(alpha_o);
else
    Fc = 0;
    fErr = [];
end
Fi = 0*Qi;
Fo = 0*Qo;

% race compliance
if B.Options.bRaceCompliancei
    Qri = -E.r*race_compliance(B.Race.Inner,-wi);
    fErr = [fErr;
            Qi - Qri];
else
    Qri = Qi;
end
if B.Options.bRaceComplianceo
    Qro = E.r*race_compliance(B.Race.Outer, wo);
    fErr = [fErr;
            Qo - Qro];
else
    Qro = Qo;
end

Xz = B.Geometry.D/2*sin(E.alpha)*x0 + db.*sin(alpha_i);
Xr = B.Geometry.D/2*cos(E.alpha)*x0 + db.*cos(alpha_i);

%% Race loads
Fri = Qi.*cos(alpha_i) + Fi.*sin(alpha_i);
Fzi = Qi.*sin(alpha_i) - Fi.*cos(alpha_i);
Wi = [sum(Fri.*cos(PSI));
    sum(Fri.*sin(PSI));
    sum(Fzi);
    sum( B.Geometry.rRacei.*Fzi.*sin(PSI) - Z.*Fri.*sin(PSI));
    sum(-B.Geometry.rRacei.*Fzi.*cos(PSI) + Z.*Fri.*cos(PSI));
    0*x0];

Fro = Qo.*cos(alpha_o) + Fo.*sin(alpha_o);
Fzo = Qo.*sin(alpha_o) - Fo.*cos(alpha_o);
Wo =-[sum(Fro.*cos(PSI));
      sum(Fro.*sin(PSI));
      sum(Fzo);
      sum( B.Geometry.rRaceo.*Fzo.*sin(PSI) - Z.*Fro.*sin(PSI));
      sum(-B.Geometry.rRaceo.*Fzo.*cos(PSI) + Z.*Fro.*cos(PSI));
      0*x0];

%forces
F.Fi = Wi;
F.Fo = Wo;
F.FInt = fErr;

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