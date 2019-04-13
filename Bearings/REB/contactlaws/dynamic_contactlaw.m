function [Qi,vr,wi,wo,K] = dynamic_contactlaw(Contact,Race,Options,Fc,dn)

[wi,wo,vr] = ball_newton(Contact,Race,Options,Fc,dn);
[Qi,Ki] = hertz_contactlaw(Contact.Inner.K,Contact.n,dn-vr-wi,Contact.tol);
[~,Ko] = hertz_contactlaw(Contact.Outer.K,Contact.n,vr-wo,Contact.tol);

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
K = 1./(1./Ki + 1./Ko + 1./Kri + 1./Kro);

function [wi,wo,vr,iter] = ball_newton(Contact,Race,Options,Fc,dn)
vr = dn / (1 + B.Contact.lambda);

%find elements out of contact
db_crit = (Fc/Contact.Outer.K).^(1/Contact.n);
iOuterOnly = dn<db_crit;
vr(iOuterOnly) = db_crit(iOuterOnly);

sz = size(dn);
dn = permute(dn(:),[2 3 1]);
vr = permute(vr(:),[2 3 1]);
Fc = permute(Fc(:),[2 3 1]);
x = [0*dn; 0*dn; vr];
dQ = x + Inf;
iter = 0;

iFlex = [];
if Options.bRaceCompliancei
    iFlex(end+1) = 1;
end
if Options.bRaceComplianceo
    iFlex(end+1) = 2;
end
if Options.bCentrifugal
    iFlex(end+1) = 3;
end

Qri = dn + NaN; Kri = dn + NaN;
Qro = dn + NaN; Kro = dn + NaN;

while any(abs(dQ(:))>1E-8)
    wi = x(1,:,:);
    wo = x(2,:,:);
    vr = x(3,:,:);
    
    if Options.bCentrifugal
        [Qi,Ki] = hertz_contactlaw(Contact.Inner.K,Contact.n,dn-vr-wi,Contact.tol);
        [Qo,Ko] = hertz_contactlaw(Contact.Outer.K,Contact.n,vr-wo,Contact.tol);
    else
        [Qi,Ki] = hertz_contactlaw(Contact.K,Contact.n,dn-wi-wo,Contact.tol);
        Qo = Qi; Ko = Ki;
    end
    
    if Options.bRaceCompliancei
        [Qri,Kri] = race_compliance_loads(Race.Inner,-wi); Qri = -Qri;
    end
    if Options.bRaceComplianceo
        [Qro,Kro] = race_compliance_loads(Race.Outer,wo);
    end
    
    Kmat = [Ki+Kri  Ki     Ki;
             Ko   Ko+Kro  -Ko;
             Ki    -Ko    Ki+Ko];

    dQ = [Qi-Qri;
          Qo-Qro;
          Qi-Qo+Fc];  

    x(iFlex,:,:) = x(iFlex,:,:) + mtimesx(minvx(Kmat(iFlex,iFlex,:)),dQ(iFlex,:,:));
    iter = iter + 1;
    if iter > 100
        break
    end
end

wi = reshape(x(1,:,:),sz);
wo = reshape(x(2,:,:),sz);
vr = reshape(x(3,:,:),sz);