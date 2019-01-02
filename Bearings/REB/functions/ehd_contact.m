function [Qi,Ki,Ci] = ehd_contact(I,us,dbi,dbidot)

n   = I.EHD.n0 + I.EHD.Cn.*(1-exp(-I.EHD.kn*us.^I.EHD.en));
db0 = I.EHD.Cd*us.^I.EHD.ed;
[Qi,Ki] = hertz_contact(I.K,n,dbi-db0);

Ci = 1;
Qi = Qi + Ci*dbidot;