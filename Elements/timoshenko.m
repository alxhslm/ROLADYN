function [Ke,Me,Ge] = timoshenko(material,ri,ro,le)
if isinf(material.E)
    E = 0;
else
    E = material.E;
end
rho = material.rho;

%shear modulus
v = material.v;
Ge = E/2/(1+v);

%derived geometric params
Ae = pi*(ro^2 - ri^2);
Ie = pi/4*(ro^4 - ri^4);
mu = ri/ro;

%shear deformation params
k = (6*(1+v)*(1+mu^2)^2) / ((7+6*v)*(1+mu^2)^2 + (20+12*v)*mu^2);
phi = 12*E*Ie/(k*Ge*Ae*le^2 + eps);

%stiffness
K = stiff_coeffs(le,phi)*E*Ie;
Ke = blkdiag(K,K);

%mass
[M_lin,M_rot] = mass_coeffs(le,phi);
M = rho*(M_lin*Ae + M_rot*Ie);
Me = blkdiag(M,M);

%gyro
G = gyro_coeffs(le,phi)*rho*Ie;
Ge = antiblkdiag(G,-G);
    
function k = stiff_coeffs(le,phi)
k  = [ 12      6*le           -12       6*le;
      6*le   (4+phi)*le^2   -6*le     (2-phi)*le^2;
      -12     -6*le           12       -6*le;    
      6*le  (2-phi)*le^2    -6*le     (4+phi)*le^2]/(1+phi)/le^3;
  
function [m_lin,m_rot] = mass_coeffs(le,phi)
m1 = 312 + 588*phi + 280*phi^2;
m2 = (44 + 77*phi + 35*phi^2)*le;
m3 = 108 + 252*phi + 140*phi^2;
m4 = -(26 + 63*phi + 35*phi^2)*le;
m5 = (8+14*phi+7*phi^2)*le^2;
m6 = -(6+14*phi+7*phi^2)*le^2;
m7 = 36;
m8 = (3-15*phi)*le;
m9 = (4+5*phi+10*phi^2)*le^2;
m10 = (-1-5*phi+5*phi^2)*le^2;

m_lin  = [m1  m2  m3  m4;
          m2  m5 -m4  m6;
          m3 -m4  m1 -m2;
          m4  m6 -m2  m5]*le/840/(1+phi)^2;

m_rot  = [m7  m8  -m7  m8;
          m8  m9  -m8  m10;
         -m7 -m8   m7 -m8;
          m8  m10 -m8  m9]/30/le/(1+phi)^2;
  
function g = gyro_coeffs(le,phi)
g1 = 36;
g2 = (3-15*phi)*le;
g3 = (4+5*phi+10*phi^2)*le^2;
g4 = (-1-5*phi+5*phi^2)*le^2;

g  = [g1   g2    -g1   g2;
     g2    g3    -g2   g4;
     -g1  -g2     g1  -g2;
     g2    g4    -g2   g3]/15/le/(1+phi)^2;
