function [Ke,Me,Ge] = euler_bernoulli(material,ri,ro,le)
if isinf(material.E)
    E = 0;
else
    E = material.E;
end
rho = material.rho;

%derived geometric params
Ae = pi*(ro^2 - ri^2);
Ie = pi/4*(ro^4 - ri^4);

%stiffness
K =  [ 12    6*le   -12     6*le;
      6*le   4*le^2 -6*le   2*le^2;
      -12   -6*le    12    -6*le;
      6*le   2*le^2 -6*le   4*le^2]*E*Ie/le^3;
	  
Ke = blkdiag(K,K);
   
%mass
M_lin = [ 156   22*le     54   -13*le;
         22*le   4*le^2  13*le  -3*le^2;  
          54    13*le     156  -22*le;
        -13*le  -3*le^2 -22*le   4*le^2]*le/420;

M_rot =  [ 36       3*le    -36      3*le; 
          3*le     4*le^2   -3*le   -le^2;
          -36     -3*le      36     -3*le; 
          3*le     -le^2    -3*le   4*le^2]/le/30; 
	  
M = rho*(M_lin*Ae + M_rot*Ie);
Me = blkdiag(M,M);
 
%gyro
G  = [36   3*le    -36   3*le
     3*le  4*le^2 -3*le  -le^2;
     -36  -3*le     36  -3*le
     3*le -le^2   -3*le  4*le^2]*rho*Ie/(15*le);
 
Ge = antiblkdiag(G,-G);