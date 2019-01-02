function [F,V,S] =  SFD_sanandres(D,States)
q    = States.qi - States.qo;
qdot = States.qidot - States.qodot;

NPts = size(q,2);

%the angle is specified w.r.t to the X axis
theta = linspace(-pi,pi,D.Nt);
z = linspace(-D.L/2,D.L/2,D.Nz);
[T,Z] = meshgrid(theta,z);

%repeat for every timestep
T = repmat(T,1,1,NPts);
Z = repmat(Z,1,1,NPts);

q(1,:) = q(1,:) + eps;
Q = repmat(permute(q,[3 4 2 1]),D.Nz,D.Nt);
qdot = repmat(permute(qdot,[3 4 2 1]),D.Nz,D.Nt);

ex = Q(:,:,:,1) + Z .* Q(:,:,:,4);
ey = Q(:,:,:,2) - Z .* Q(:,:,:,3);

e = sqrt(ex.^2 + ey.^2);
er = cat(4,ex./e,ey./e);
et = cat(4,-er(:,:,:,2),er(:,:,:,1));

edot  = cat(4, qdot(:,:,:,1) + Z .* qdot(:,:,:,4), ...
               qdot(:,:,:,2) - Z .* qdot(:,:,:,3));
          
vt = sum(et.*edot,4);
omega = vt./e;

T = T - atan2(ey,ex) - pi;
p = sanandres(e,omega,T,Z,D) + D.pSupply;
p = max(p,D.pCavitation);
T = T + atan2(ey,ex) + pi;

%our theta trick means integrating the forces is super easy
W = [permute(trapz(z,trapz(theta, D.R*p.*cos(T),2),1),[1 3 2]);
     permute(trapz(z,trapz(theta, D.R*p.*sin(T),2),1),[1 3 2]);
     0*wons;
     permute(trapz(z,trapz(theta,-D.R*p.*sin(T).*Z,2),1),[1 3 2]);
     permute(trapz(z,trapz(theta, D.R*p.*cos(T).*Z,2),1),[1 3 2]);
     0*wons];
 
%forces
F.F = W;
F.FInt = zeros(0,NPts);

%channels
V.p = p;
V.T = T;
V.Z = Z;
V.vr = vr;
V.vt = vt;
V.ar = ar;
V.at = at;

%damping
S = struct();

function p = sanandres(e,omega,theta,z,P)

sigma = e/P.c;
delta = P.c/P.R;
Re = omega.*P.c^2*P.rho/P.mu;

zeta = (1-1i)*sqrt(Re/2);

xi = z./P.R;

Z = 1 - cosh(xi)/cosh(0.5*P.L/P.R);

k = zeta.* sinh(zeta) ./ (zeta.*h.*sinh(zeta) + 2 - 2*cosh(zeta) + eps);

p_bar = sigma .* Re .* Z .* real(k.*exp(1i*theta));
p = p_bar .* omega * P.mu / delta^2;