function [Qi,Ki] = hertz_contact(Ki,n,dbi,tol)
if nargin < 4
    tol = 0;
end
if nargout > 1
    [dbi,d2] = maxSmooth(dbi,0,tol);
    [Qi,d1] = sgn_power(dbi,n);
    Qi = Qi * Ki;
%     [Qi,d2] = maxSmooth(Qi,0,tol);
    Ki = Ki .* d1 .* d2;
else
    dbi = maxSmooth(dbi,0,tol);
    Qi = Ki*sgn_power(dbi,n);
%     Qi = maxSmooth(Qi,0,tol);
end