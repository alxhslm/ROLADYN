function [F,V,S] = piezo_goldfarb(B,States)

xPz    = States.qi(1,:)    - States.qo(1,:);
xPzdot = States.qidot(1,:) - States.qodot(1,:);
xPzddot = States.qiddot(1,:) - States.qoddot(1,:);

z    = States.xInt;
zdot = States.xIntdot;

Fm    = B.Mech.k*xPz    + B.Mech.c*xPzdot;
Fmdot = B.Mech.k*xPzdot + B.Mech.c*xPzddot;

q    = B.T*xPz    + B.C/B.T*Fm;
qdot = B.T*xPzdot + B.C/B.T*Fmdot;

if nargout < 2
    zdot_model = bouc_wen_ode(B.Bouc_wen,qdot,z);
else
    [zdot_model,dzdot_dz,dzdot_dqdot] = bouc_wen_ode(B.Bouc_wen,qdot,z);
end

F.F = B.Mech.k*xPz + B.Mech.c*xPzdot - B.T*z/B.Cm;
F.FInt = (zdot - zdot_model);

V.q = q;
V.qdot = qdot;
V.Vpz = Fm/B.T;
V.Vh  = z/B.Cm;

if nargout > 1
    S.Kqq =  B.Mech.k;
    S.Kqx = -B.T/B.Cm;
    S.Kxq =  0;
    S.Kxx = -dzdot_dz;
    
    S.K = -B.kS;
    
    S.Cqq = B.Mech.c;
    S.Cqx =  0;
    S.Cxq = -dzdot_dqdot * (B.T + B.C/B.T*B.Mech.k);
    S.Cxx = 1;
    S.C = -B.Mech.c;
    
    S.Mqq = 0;
    S.Mqx = 0;
    S.Mxq = -dzdot_dqdot * (B.C/B.T*B.Mech.c);
    S.Mxx = 0;
    S.M   = 0;
end

function [zdot,dzdot_dz,dzdot_dqdot] = bouc_wen_ode(p,qdot,z)
z = z + eps;
zdot = qdot  - p.beta*abs(qdot).*(abs(z).^(p.n-1)).*z - p.gamma*qdot.*(abs(z).^p.n);
dzdot_dqdot = 1 - p.beta*sign(qdot).*(abs(z).^(p.n-1)).*z - p.gamma.*(abs(z).^p.n);
dzdot_dz = - p.beta*abs(qdot).*((abs(z).^(p.n-1)) + (p.n-1)*(abs(z).^(p.n-2)).*sign(z).*z) - p.n*p.gamma.*(abs(z).^(p.n-1)).*sign(z);