function [F,V,S] = piezo_goldfarb(B,States)

xPz     = States.qi     - States.qo;
xPzdot  = States.qidot  - States.qodot;

Vpz     = States.xInt(1,:);
Vpzdot  = States.xIntdot(1,:);

z    = States.xInt(2,:);
zdot = States.xIntdot(2,:);

Ve = B.G*States.u;

q    = B.T*xPz    + B.C*Vpz;
qdot = B.T*xPzdot + B.C*Vpzdot;

if nargout < 2
    zdot_model = bouc_wen_ode(B.Bouc_wen,qdot,z);
else
    [zdot_model,dzdot_dz,dzdot_dqdot] = bouc_wen_ode(B.Bouc_wen,qdot,z);
end

Vh = z/B.Cm;

F.F = B.Mech.k*xPz + B.Mech.c*xPzdot - B.T*Vpz;

F.FInt = [Vpz + Vh - Ve;
          (zdot - zdot_model)];

V.q    = q;
V.qdot = qdot;

V.z    = z;
V.zdot = zdot;

V.Vpz  = Vpz;
V.Vh   = Vh;
V.Ve   = Ve;

V.Fpz  = B.T*Vpz;

if nargout > 1
    wons = permute(0*xPz,[1 3 2])+1;
    S.Kqq =  B.Mech.k*wons;
    S.Kqx = [-B.T*wons 0*wons];
    S.Kxq = [0*wons;
             0*wons];
    S.Kxx =  [wons       1/B.Cm*wons;
              0*wons   -permute(dzdot_dz,[1 3 2])];
    S.Kqu = 0*wons;
    S.Kxu = [-B.T*B.G*wons;
             0*wons];
    
    S.Cqq = B.Mech.c*wons;
    S.Cqx = [0*wons 0*wons];
    S.Cxq = [0*wons;
             -permute(dzdot_dqdot .* B.T,[1 3 2])];
    S.Cxx = [0*wons    0*wons;
             -permute(dzdot_dqdot .* B.C,[1 3 2]) wons];
    S.Cqu = 0*wons;
    S.Cxu = [0*wons;
             0*wons];
end

function [zdot,dzdot_dz,dzdot_dqdot] = bouc_wen_ode(p,qdot,z)
zdot = qdot  - p.beta*abs(qdot).*(abs(z).^(p.n-1)).*z - p.gamma*qdot.*(abs(z).^p.n);
dzdot_dqdot = 1 - p.beta*sign(qdot).*(abs(z).^(p.n-1)).*z - p.gamma.*(abs(z).^p.n);
dzdot_dz = - p.beta*abs(qdot).*((abs(z).^(p.n-1)) + (p.n-1)*(abs(z).^(p.n-2)).*sign(z).*z) - p.n*p.gamma*qdot.*(abs(z).^(p.n-1)).*sign(z);