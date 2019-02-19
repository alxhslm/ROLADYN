function [F,V,S] =  SFD_tichy(D,States)
q    = States.qi - States.qo;
qdot = States.qidot - States.qodot;
qddot = States.qiddot - States.qoddot;
NPts = size(q,2);

%the angle is specified w.r.t to the X axis
theta = linspace(-pi,pi,D.Nt);
z = linspace(-D.L/2,D.L/2,D.Nz);
[T,Z] = meshgrid(theta,z);

%repeat for every timestep
T = repmat(T,1,1,NPts);
Z = repmat(Z,1,1,NPts);

wons = q(1,:)*0 + 1;
q(1,:) = q(1,:) + eps;
Q = repmat(permute(q,[3 4 2 1]),D.Nz,D.Nt);
qdot = repmat(permute(qdot,[3 4 2 1]),D.Nz,D.Nt);
qddot = repmat(permute(qddot,[3 4 2 1]),D.Nz,D.Nt);

if D.bOutOfPlane
    Zoffset = Z;
else 
    Zoffset = 0*Z;
end

ex = Q(:,:,:,1) + Zoffset .* Q(:,:,:,5);
ey = Q(:,:,:,2) - Zoffset .* Q(:,:,:,4);

e = sqrt(ex.^2 + ey.^2);
er = cat(4,ex./e,ey./e);
et = cat(4,-er(:,:,:,2),er(:,:,:,1));

edot  = cat(4, qdot(:,:,:,1) + Zoffset .* qdot(:,:,:,5), ...
               qdot(:,:,:,2) - Zoffset .* qdot(:,:,:,4));
          
eddot = cat(4, qddot(:,:,:,1) + Zoffset .* qddot(:,:,:,5), ...
               qddot(:,:,:,2) - Zoffset .* qddot(:,:,:,4));

vr = sum(er.*edot,4);
vt = sum(et.*edot,4);
ar = sum(er.*eddot,4);
at = sum(et.*eddot,4);

%now make the angle relative the shaft centre at each axial location
T = T - atan2(ey,ex) - pi;
p = reynolds(e,vr,vt,ar,at,T,Z,D) + D.pSupply;
p = max(p,D.pCavitation);

%put the angle back to fixed coords
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
 
function p = reynolds(e,vr,vt,ar,at,theta,z,P)
h = P.c + e.*cos(theta);

p = P.mu  ./ (h.^3) .*  12  .* (vr.*cos(theta) + vt.*sin(theta)) + ...
    P.rho ./  h     .*  1.2 .* (ar.*cos(theta) + at.*sin(theta));
p = p.*(z.^2 - (P.L/2).^2)/2;