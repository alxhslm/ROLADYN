function [Ke,Me,Ge,Kc] = annular(material,t,th,r,beta,dR)
if isinf(material.E)
    E = 0;
else
    E = material.E;
end
v = material.v;
tensor = E/(1-v^2) * [ 1     v        0;
                       v     1        0;
                       0     0     0.5*(1-v)];
rho = material.rho;

NInt = [10 12];
uth = linspace(0,beta,NInt(1)+1);
dth = diff(uth);
uth = 0.5*(uth(1:end-1)+uth(2:end));

ur = linspace(0,dR,NInt(2)+1);
dr = diff(ur);
ur = 0.5*(ur(1:end-1)+ur(2:end));

[uth,ur] = meshgrid(uth,ur);
uth = permute(uth(:),[2 3 1]); ur = permute(ur(:),[2 3 1]);

[dth,dr] = meshgrid(dth,dr);
dth = permute(dth(:),[2 3 1]); dr = permute(dr(:),[2 3 1]);

r = r + ur;
th = th + uth;

%extract displacements/derivatives from shape fun
C = quad_rect_shapefun('coeff',beta,dR);
[P,Deriv,Deriv2] = quad_rect_shapefun({'fun','jac','hess'},uth,ur);
invC = inv(C);

dA = (r .* dth .* dr);

%first handle the mass matrix
%w = P*coeffs
%dw = Q*coeffs
N =  mtimesx(P,invC);
Mlat = rho*t*sum(mtimesx(mtransposex(N),N).*dA,3);

dw_theta = Deriv(1,:,:); Bt  = mtimesx(dw_theta./r,invC);
dw_dr    = Deriv(2,:,:); Br  = mtimesx(dw_dr,invC);

BB = mtimesx(mtransposex(Bt),Bt) + mtimesx(mtransposex(Br),Br);
Mrot =  rho*t^3/12*sum(BB.*dA,3);
Me = Mlat + Mrot;

%gyroscopic matrix
NB = mtimesx(mtransposex(Bt),N) - mtimesx(mtransposex(N),Bt);
Ge = -rho*t*sum(r.*NB.*dA,3);

zer = 0*th;
won = zer + 1;
Strain = [zer     1./r  1/r.^2 zer  zer;
          zer     zer   zer    won  zer;
         -1./r.^2 zer   zer    zer 1./r];

%strain matrix
Qe = mtimesx(Strain,[Deriv; Deriv2]);

D = (t^3)/12 * tensor;
                        
B  = mtimesx(Qe,invC);
V  = mtimesx(mtransposex(B),mtimesx(D,B));

Ke = sum(V.*dA,3);

%centrifugal stiffnening
Kc   =  Mrot;