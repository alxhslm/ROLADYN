function D = mass2dim(D)
dims = {'Ri';
        'Ro';
        't'};

Geom.Ri = D.Ring.R(1);
Geom.Ro = D.Ring.R(2);
Geom.t  = D.Ring.t;
    
bKnownX = false(3,1);
Xknown = zeros(3,1);        
for i = 1:size(dims,1)
    if isfield(Geom,dims{i})
        field = Geom.(dims{i});
        if ~isnan(field)
            bKnownX(i) = true;
            Xknown(i) = field;
        end
    end
end

inta = {'m';
        'Id';
        'Ip'};
    
bKnownF = false(3,1);
Fknown = zeros(3,1);        
for i = 1:size(inta,1)
    if isfield(D.Inertia,inta{i})
        field = D.Inertia.(inta{i});
        if ~isnan(field)
            bKnownF(i) = true;
            Fknown(i) = field;
        end
    end
end

if (sum(bKnownF) + sum(bKnownX)) < 3
    error('Not enough information given for disc "%s"',D.Name);
end

if (sum(bKnownF) + sum(bKnownX)) > 3
    error('Too much information given for disc "%s"',D.Name);
end

Xscale = [1E-3;1E-2;1E-2];
fscale = [1;1E-4;1E-4];

Xscale = Xscale(~bKnownX);
fscale = fscale(bKnownF);

if ~all(bKnownX)
    X0 = rand(sum(~bKnownX),1);
    constr_fun = @(X)confun(X,Xscale,fscale,Xknown,Fknown,bKnownX,bKnownF,D.Material.rho);
    options.print_level = 0;
    iter = 1;
    while iter < 10
        [X,info] = fipopt([],X0,constr_fun,options);
        X0 = X;
        iter = iter + 1;
    end
    if info.status
         error('Unable to find disc parameters')
    end
else
    X = [];
end

Xsol = Xknown;
Xsol(~bKnownX) = X.*Xscale;

[ri,ro,t] = X2props(Xsol);
[m,Id,Ip] = disc_properties(D.Material.rho,ri,ro,t);
f = [m;Id;Ip];

for i = 1:size(dims,1)
    if ~bKnownX(i)
        Geom.(dims{i}) = Xsol(i);
    end
end

D.Ring.R = [Geom.Ri Geom.Ro];
D.Ring.t = Geom.t;

for i = 1:size(inta,1)
    if ~bKnownF(i)
        D.Inertia.(inta{i}) = f(i);
    end
end

function err = confun(X,Xscale,fscale,Xknown,Fknown,bKnownX,bKnownF,rho)
Xtest = Xknown;
Xtest(~bKnownX) = X.*Xscale;
[ri,ro,t] = X2props(Xtest);
[m,Id,Ip] = disc_properties(rho,ri,ro,t);
f = [m;Id;Ip];
err = f - Fknown;
err = err(bKnownF) ./ fscale;

function [ri,ro,t] = X2props(X)
ri   = X(1);
ro   = X(2);
t    = X(3);