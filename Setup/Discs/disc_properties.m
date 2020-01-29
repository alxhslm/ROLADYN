function [m,Id,Ip] = disc_properties(rho,ri,ro,t)
ri2 = ri^2 * sign(ri);
ro2 = ro^2 * sign(ro);
ri4 = ri^2 * ri2;
ro4 = ro^2 * ro2;

m  = rho*pi*t*(ro2 - ri2);
if nargout > 1
    Id = rho*pi*t/12*(3*(ro4 - ri4) + t^2*(ro2-ri2));
end
if nargout > 2
    Ip = rho*pi*t/2*(ro4 - ri4);
end