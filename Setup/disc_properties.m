function [m,Id,Ip] = disc_properties(rho,ri,ro,t)
m  = rho*pi*t*(ro^2 - ri^2);
Id = rho*pi*t/12*(3*(ro^4 - ri^4) + t^2*(ro^2-ri^2));
Ip = rho*pi*t/2*(ro^4 - ri^4);