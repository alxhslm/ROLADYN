function [F,V,S] = REB_ehd(B,States)

qi = States.qi;
qo = States.qo;
qidot = States.qidot;
qodot = States.qodot;

xInt = States.xInt;
xdotInt = States.xdotInt;
xddotInt = States.xddotInt;

Oi = States.Oi;
Oo = States.Oo;
Ai = States.Ai;
Ao = States.Ao;

E = B.Elements;

if B.Options.bSjovall
    E = B.Sjovall;
else
    E = B.Balls;
end

%create some vectors of the right size
wons = (E.psi*0+1);
x0 = (qi(1,:)*0 + 1);

%if we're not including VC effects, set Ocage to 0 to "fix" the cage
if ~B.Options.bVC
    Acage = 0*Ai;
    Ocage = 0*Oi;
else
    Acage = B.Kinematics.rCagei * Ai + B.Kinematics.rCageo * Ao;
    Ocage = B.Kinematics.rCagei * Oi + B.Kinematics.rCageo * Oo;
end

vz = xInt(1:E.N,:);
vr = xInt(E.N+(1:E.N),:);

vzdot = xdotInt(1:E.N,:);
vrdot = xdotInt(E.N+(1:E.N),:);

vzddot = xddotInt(1:E.N,:);
vrddot = xddotInt(E.N+(1:E.N),:);

z = B.Geometry.zRacei*sign(E.z/(E.z(1)+eps));
PSI = E.psi*x0 + wons*Acage;
Z = z*x0;

dzi  = wons*qi(3,:) + B.Geometry.rRacei*(sin(PSI).*(wons*qi(4,:)) - cos(PSI).*(wons*qi(5,:)));
dri  = cos(PSI).*(wons*qi(1,:)) + sin(PSI).*(wons*qi(2,:)) - Z.*(sin(PSI).*(wons*qi(4,:)) - cos(PSI).*(wons*qi(5,:)));

dzo  = wons*qo(3,:) + B.Geometry.rRaceo*(sin(PSI).*(wons*qo(4,:)) - cos(PSI).*(wons*qo(5,:)));
dro  = cos(PSI).*(wons*qo(1,:)) + sin(PSI).*(wons*qo(2,:)) - Z.*(sin(PSI).*(wons*qo(4,:)) - cos(PSI).*(wons*qo(5,:)));

dzdoti  = wons*qidot(3,:) + B.Geometry.rRacei*(sin(PSI).*(wons*qidot(4,:)) - cos(PSI).*(wons*qidot(5,:))) + ...
         Ocage*B.Geometry.rRacei*(cos(PSI).*(wons*qi(4,:)) + sin(PSI).*(wons*qi(5,:)));

drdoti  = cos(PSI).*(wons*qidot(1,:)) + sin(PSI).*(wons*qidot(2,:)) - Z.*(sin(PSI).*(wons*qidot(4,:)) - cos(PSI).*(wons*qidot(5,:))) + ... 
        Ocage*(-sin(PSI).*(wons*qi(1,:)) + cos(PSI).*(wons*qi(2,:)) - Z.*(cos(PSI).*(wons*qi(4,:)) + sin(PSI).*(wons*qi(5,:))));
    
dzdoto  = wons*qodot(3,:) + B.Geometry.rRaceo*(sin(PSI).*(wons*qodot(4,:)) - cos(PSI).*(wons*qodot(5,:))) + ...
         Ocage*B.Geometry.rRacei*(cos(PSI).*(wons*qo(4,:)) + sin(PSI).*(wons*qo(5,:)));     
    
drdoto  = cos(PSI).*(wons*qodot(1,:)) + sin(PSI).*(wons*qodot(2,:)) - Z.*(sin(PSI).*(wons*qodot(4,:)) - cos(PSI).*(wons*qodot(5,:))) + ... 
        Ocage*(-sin(PSI).*(wons*qo(1,:)) + cos(PSI).*(wons*qo(2,:)) - Z.*(cos(PSI).*(wons*qo(4,:)) + sin(PSI).*(wons*qo(5,:))));

Az = B.Geometry.A0*sin(E.alpha)*x0 + dzi - dzo - B.Geometry.cz;
Ar = B.Geometry.A0*cos(E.alpha)*x0 + dri - dro - B.Geometry.cr;

%% Contact forces
Xz = vz - dzo;
Xr = vr - dro;
[Ai,Ao,alpha_i,alpha_o] = race_geometry(Xz,Xr,Az,Ar);
dbi = Ai - (B.Geometry.RRacei-B.Geometry.D/2);
dbo = Ao - (B.Geometry.RRaceo-B.Geometry.D/2); 
dbi_dot = (dzdoti-vzdot).*sin(alpha_i) + (drdoti - vrdot).*cos(alpha_i);
dbo_dot =-(dzdoto-vzdot).*sin(alpha_i) - (drdoto - vrdot).*cos(alpha_o);
[usi,uso] = ehd_sum_speeds(B,alpha_i,alpha_o,wons*Oi,wons*Oo);
Qi = E.r*ehd_contactlaw(B.Contact.Inner,usi,dbi,dbi_dot);
Qo = E.r*ehd_contactlaw(B.Contact.Outer,uso,dbo,dbo_dot);
[Fc,Fi,Fo,Mg] = dynamic_ball_loads(B,alpha_i,alpha_o,wons*Oi,wons*Oo);
Fc = E.r*Fc; Fi = E.r*Fi; Fo = E.r*Fo; Mg = E.r*Mg;
fErr = [Qi.*sin(alpha_i) - Fi.*cos(alpha_i) - Qo.*sin(alpha_o) + Fo.*cos(alpha_o)      - B.Dynamics.mb * vzddot;
        Qi.*cos(alpha_i) + Fi.*sin(alpha_i) - Qo.*cos(alpha_o) - Fo.*sin(alpha_o) + Fc - B.Dynamics.mb * vrddot];

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
Wo = [sum(Fro.*cos(PSI));
    sum(Fro.*sin(PSI));
    sum(Fzo);
    sum( B.Geometry.rRaceo.*Fzo.*sin(PSI) - Z.*Fro.*sin(PSI));
    sum(-B.Geometry.rRaceo.*Fzo.*cos(PSI) + Z.*Fro.*cos(PSI));
    0*x0];


%forces
F.F = Wi;
F.FInt = fErr;

%contact loads
V.Qi = Qi; V.Qo = Qo; 
V.Fi = Fi; V.Fo = Fo;

%geometry
V.alpha_i = alpha_i; V.alpha_o = alpha_o;
V.dbi = dbi; V.dbo = dbo;
V.Ai = Ai;  V.Ao = Ao;
V.Xr = Xr;  
V.Xz = Xz;

%dynamic loads
V.Fc = Fc;   
V.Mg = Mg;

%stiffnesses
S = struct([]);