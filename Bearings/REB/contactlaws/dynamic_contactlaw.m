function [Qi,Qo,u,Ki,Ko] = dynamic_contactlaw(Contact,Fc,r,tol)
if nargin < 4
    tol = 0;
end

lambda = (Contact.Outer.K / Contact.Inner.K)^(1/Contact.n);
u = r / (1 + lambda);
dc = (Fc/Contact.Outer.K).^(1/Contact.n);
iO = r<dc;
u(iO) = dc(iO);
err = u + Inf;
iter = 0;
while any(abs(err(:))>1E-8)
    [Qi,Ki] = hertz_contactlaw(Contact.Inner.K,Contact.n,r-u,tol);
    [Qo,Ko] = hertz_contactlaw(Contact.Outer.K,Contact.n,u,tol);
    err = Qi + Fc - Qo;
    u = u + err./(Ki + Ko);
    iter = iter + 1;
    if iter > 50
        Qi = Qi + NaN;
        Qo = Qo + NaN;
        break;
    end
end