function [Qi,Qo] = dynamic_contact(Contact,Fc,r,tol)
if nargin < 4
    tol = 0;
end

% p = interp1q(Contact.O,Contact.p,O')';
% r = r + p(3,:);
% Qi = p(1,:).*sgn_power(r,p(2,:));
% Qi = maxSmooth(Qi,0,tol);

lambda = (Contact.Outer.K / Contact.Inner.K)^(1/Contact.n);
u = r / (1 + lambda);
dc = (Fc/Contact.Outer.K).^(1/Contact.n);
iO = r<dc;
u(iO) = dc(iO);
err = u + Inf;
if isnan(err)
    Qi = NaN;
    Qo = NaN;
    return;
end
iter = 0;
while any(abs(err(:))>1E-8)
    [Qi,Ki] = hertz_contact(Contact.Inner.K,Contact.n,r-u,tol);
    [Qo,Ko] = hertz_contact(Contact.Outer.K,Contact.n,u,tol);
    err = Qi + Fc - Qo;
    u = u + err./(Ki + Ko);
    iter = iter + 1;
    if iter > 50
        Qi = Qi + NaN;
        Qo = Qo + NaN;
        break;
    end
end