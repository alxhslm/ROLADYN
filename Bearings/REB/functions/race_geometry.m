function [li,lo,alpha_i,alpha_o] = race_geometry(Xz,Xr,Az,Ar)
%compute contact angle and deflections
li = sqrt((Az-Xz).^2  + (Ar-Xr).^2);
lo = sqrt(Xz.^2  + Xr.^2);

alpha_i = atan2(Az-Xz,Ar-Xr);
alpha_o = atan2(Xz,Xr);
