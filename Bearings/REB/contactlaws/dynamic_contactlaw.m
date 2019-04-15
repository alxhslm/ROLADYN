function [Qi,Qo,vr,wi,wo,K] = dynamic_contactlaw(Contact,Race,Options,Fc,dn)

[wi,wo,vr] = ball_newton(Contact,Race,Options,Fc,dn);
[~,~,Qi,Qo,Kri,Kro,Ki,Ko] = get_forces(Contact,Race,Options,dn,wi,wo,vr);
K = 1./(1./Ki + 1./Ko + 1./Kri + 1./Kro);

function [wi,wo,vr,iter] = ball_newton(Contact,Race,Options,Fc,dn)
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
vr = max(dn,0) / (1 + Contact.lambda);

%find elements out of contact
dn_crit = (Fc/Contact.Outer.K).^(1/Contact.n);
iOuterOnly = dn<dn_crit;
wi(iOuterOnly) = 0;
wo(iOuterOnly) = (Fc(iOuterOnly)/Kro);
vr(iOuterOnly) = dn_crit(iOuterOnly);

%make everything 1D
sz = size(dn);
dn = permute(dn(:),[2 3 1]);
Fc = permute(Fc(:),[2 3 1]);

wi = permute(wi(:),[2 3 1]);
wo = permute(wo(:),[2 3 1]);
vr = permute(vr(:),[2 3 1]);

iFlex = [];
if Options.bRaceCompliancei
    iFlex(end+1) = 1;
end
if Options.bRaceComplianceo
    iFlex(end+1) = 2;
end
iFlex(end+1) = 3;

x = [wi; wo; vr];
dQ = x + Inf;
iter = 0;
Z = 0*dn;

while any(abs(dQ(:))>1E-8)
    wi = x(1,:,:);
    wo = x(2,:,:);
    vr = x(3,:,:);
    
    [Qri,Qro,Qi,Qo,Kri,Kro,Ki,Ko] = get_forces(Contact,Race,Options,dn,wi,wo,vr);
         
    Kmat = [Ki+Kri  Z     Ki;
             Z   Ko+Kro  -Ko;
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

function [Qri,Qro,Qi,Qo,Kri,Kro,Ki,Ko] = get_forces(Contact,Race,Options,dn,wi,wo,vr)
[Qi,Ki] = hertz_contactlaw(Contact.Inner.K,Contact.n,dn-vr-wi,Contact.tol);
[Qo,Ko] = hertz_contactlaw(Contact.Outer.K,Contact.n,vr-wo,Contact.tol);

if Options.bRaceCompliancei
    [Qri,Kri] = race_compliance_loads(Race.Inner,-wi); Qri = -Qri;
else
    Qri = Qi;
    Kri = dn+Inf;
end
if Options.bRaceComplianceo
    [Qro,Kro] = race_compliance_loads(Race.Outer,wo);
else
    Qro = Qo;
    Kro = dn+Inf;
end