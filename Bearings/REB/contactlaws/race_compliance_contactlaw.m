function [Qi,Qo,wi,wo,Ktoti,Ktoto] = race_compliance_contactlaw(Contact,Race,Options,dn)

[wi,wo] = ball_newton(Contact,Race,Options,dn);

[Qi,Kh] = hertz_contactlaw(Contact.K,Contact.n,dn-wo-wi,Contact.tol);
Qo = Qi;
if Options.bRaceCompliancei
    [~,Kri] = race_compliance_loads(Race.Inner,-wi);
else
    Kri = Inf;
end
if Options.bRaceComplianceo
    [~,Kro] = race_compliance_loads(Race.Outer,wo);
else
    Kro = Inf;
end

Ktoti = 1./(1./Kh + 1./Kri + 1./Kro);
Ktoto = Ktoti;

function [wi,wo,iter] = ball_newton(Contact,Race,Options,dn)
if Options.bRaceCompliancei
    Kri = Race.Inner.K;
else
    Kri = Inf;
end
if Options.bRaceCompliancei
    Kro = Race.Outer.K;
else
    Kro = Inf;
end

%sensible initial guesses
Q0 = hertz_contactlaw(Contact.K,Contact.n,dn,Contact.tol);
wi = Q0/Kri;
wo = Q0/Kro;

%make everything 1D
sz = size(dn);
dn = permute(dn(:),[2 3 1]);
wi = permute(wi(:),[2 3 1]);
wo = permute(wo(:),[2 3 1]);

x = [wi; wo];
dQ = x + Inf;
iter = 0;

iFlex = [];
if Options.bRaceCompliancei
    iFlex(end+1) = 1;
end
if Options.bRaceComplianceo
    iFlex(end+1) = 2;
end

Qri = dn + NaN; Kri = dn + NaN;
Qro = dn + NaN; Kro = dn + NaN;

while any(abs(dQ(:))>1E-8)
    wi = x(1,:,:);
    wo = x(2,:,:);
    
    [Qi,Ki] = hertz_contactlaw(Contact.K,Contact.n,dn-wi-wo,Contact.tol);
    Qo = Qi; Ko = Ki;

    if Options.bRaceCompliancei
        [Qri,Kri] = race_compliance_loads(Race.Inner,-wi); Qri = -Qri;
    end
    if Options.bRaceComplianceo
        [Qro,Kro] = race_compliance_loads(Race.Outer,wo);
    end

    Kmat = [Ki+Kri Ki;
             Ko   Ko+Kro];

    dQ = [Qi-Qri;
          Qo-Qro];  

    x(iFlex,:,:) = x(iFlex,:,:) + mtimesx(minvx(Kmat(iFlex,iFlex,:)),dQ(iFlex,:,:));
    iter = iter + 1;
    if iter > 100
        break
    end
end

wi = reshape(x(1,:,:),sz);
wo = reshape(x(2,:,:),sz);