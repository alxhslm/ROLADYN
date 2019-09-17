function D = mass2dim(D)
fields_for_inertia = {'Geometry','Ri';
                      'Geometry','Ro';
                      'Geometry','t';
                      'Inertia','m';
                      'Inertia','Id';
                      'Inertia','Ip'};
bKnown = false(6,1);
Xknown = zeros(6,1);        
for i = 1:size(fields_for_inertia,1)
    if isfield(D.(fields_for_inertia{i,1}),fields_for_inertia{i,2})
        field = D.(fields_for_inertia{i,1}).(fields_for_inertia{i,2});
        if ~isnan(field)
            bKnown(i) = true;
            Xknown(i) = field;
        end
    end
end

Xscale = [1E-3;1E-3;1E-3;1;1E-6;1E-6];
fscale = [1;1E-6;1E-6];

Xscale = Xscale(~bKnown);

if ~all(bKnown)
    X0 = rand(sum(~bKnown),1);
    [X,~,flag] = fsolve(@(X)objfun(X,Xscale,fscale,Xknown,bKnown,D.Material.rho),X0,optimoptions('fsolve','Display','none','StepTolerance',1E-14,'FunctionTolerance',1E-14));
    if flag ~=1
        X = X + NaN;
    end
else
    X = [];
end

Xsol = Xknown;
Xsol(~bKnown) = X.*Xscale;

[D.Geometry.Ri,D.Geometry.Ro,D.Geometry.t] = X2props(Xsol);
[D.Inertia.m,D.Inertia.Id,D.Inertia.Ip] = disc_properties(D.Material.rho,D.Geometry.Ri,D.Geometry.Ro,D.Geometry.t);

    
function f = objfun(X,Xscale,fscale,Xknown,bKnown,rho)
Xtest = Xknown;
Xtest(~bKnown) = X.*Xscale;

[ri,ro,t,m0,Id0,Ip0] = X2props(Xtest);
[m,Id,Ip] = disc_properties(rho,ri,ro,t);

f = [m  -  m0;
     Id - Id0;
     Ip - Ip0];
 
 f = f ./ fscale;

function [ri,ro,t,m,Id,Ip] = X2props(X)
ri   = X(1);
ro   = X(2);
t    = X(3);
m    = X(4);
Id   = X(5);
Ip   = X(6);