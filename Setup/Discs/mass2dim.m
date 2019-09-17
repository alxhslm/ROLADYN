function D = mass2dim(D)
dims = {'Ri';
        'Ro';
        't'};

bKnownX = false(3,1);
Xknown = zeros(3,1);        
for i = 1:size(dims,1)
    if isfield(D.Geometry,dims{i})
        field = D.Geometry.(dims{i});
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

Xscale = [1E-3;1E-2;1E-2];
fscale = [1;1E-4;1E-4];

Xscale = Xscale(~bKnownX);
fscale = fscale(bKnownF);

if ~all(bKnownX)
    X0 = rand(sum(~bKnownX),1);
    [X,~,flag] = fmincon(@(X)0,X0,[],[],[],[],0*X0,X0+Inf,@(X)confun(X,Xscale,fscale,Xknown,Fknown,bKnownX,bKnownF,D.Material.rho),optimoptions('fmincon','Display','iter','StepTolerance',1E-14,'FunctionTolerance',1E-14));
    if flag ~=1
         X = X + NaN;
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
        D.Geometry.(dims{i}) = Xsol(i);
    end
end

for i = 1:size(inta,1)
    if ~bKnownF(i)
        D.Inertia.(inta{i}) = f(i);
    end
end

function [dummy,err] = confun(X,Xscale,fscale,Xknown,Fknown,bKnownX,bKnownF,rho)
Xtest = Xknown;
Xtest(~bKnownX) = X.*Xscale;
[ri,ro,t] = X2props(Xtest);
[m,Id,Ip] = disc_properties(rho,ri,ro,t);
f = [m;Id;Ip];
err = f - Fknown;
err = err(bKnownF) ./ fscale;
dummy = [];

function [ri,ro,t] = X2props(X)
ri   = X(1);
ro   = X(2);
t    = X(3);