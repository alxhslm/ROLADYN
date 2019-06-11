function D = mass2dim(D)

if D.m == 0
    D.R = [0 0.1];
    D.t = 0;
    return;
end

if ~isfield(D,'R')
    ri = 0.001;
    ro = sqrt(2*D.Ip/D.m);
    t = sqrt(12*D.Ip/D.m - 3*ro^2);

    X0 = [ri;ro;t];
    R0 = ro;
    X = fsolve(@(X)objfun(X,D,R0),X0,optimoptions('fsolve','StepTolerance',1E-14,'FunctionTolerance',1E-14));
    
    ri = X(1);
    ro = X(2);
    t =  X(3);
else
    ri = D.R(1);
    ro = D.R(2);
    t = D.m/(D.Material.rho*pi*(ro^2 - ri^2));
end

[D.m,D.Id,D.Ip] = disc_properties(D.Material.rho,ri,ro,t);
D.R = [ri ro];
D.t =  t;
    


function f = objfun(X,D,R0)
ri = X(1);
ro = X(2);
t =  X(3);

[m,Id,Ip] = disc_properties(D.Material.rho,ri,ro,t);

f(1) = m  - D.m;
f(2) = (Id - D.Id)/R0^2;
f(3) = (Ip - D.Ip)/R0^2;

