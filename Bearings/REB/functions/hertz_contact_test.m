Ki = 1;
n = 1.5;
dbi = linspace(-1,1,1001);
tol = 1E-6;
Q0 = hertz_contact(Ki,n,dbi,0);
[Q,K] = hertz_contact(Ki,n,dbi,tol);
h = 1E-8;
dQ = (hertz_contact(Ki,n,dbi+h,tol)-Q)./h;
figure
subplot(2,1,1)
plot(dbi,Q0,dbi,Q)
subplot(2,1,2)
plot(dbi,dQ)
hold on
plot(dbi,K)