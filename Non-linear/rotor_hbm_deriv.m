function Dw = rotor_hbm_deriv(dfhbm_dw,dfalg_dw,hbm,problem,w0)

c0 = hbm_nonlinear3d('func',hbm,problem,w0,x,u);
h = 1E-10;
w = w0;
for i = 1:2
    w(i) = w(i) + h;
    c = hbm_nonlinear3d('func',hbm,problem,w,x,u);
    Dw{i} = (c-c0)./h;
    w = w0;
end