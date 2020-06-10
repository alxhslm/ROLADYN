function [wi,wo] = race_compliance_contactlaw_harris(B,Ar,Az)

[wi,wo] = ball_newton(B,Ar,Az);

function [wi,wo,iter] = ball_newton(B,Ar,Az)
if B.Options.bRaceCompliancei
    Kri = B.Race.Inner.K;
else
    Kri = Inf;
end
if B.Options.bRaceCompliancei
    Kro = B.Race.Outer.K;
else
    Kro = Inf;
end

%sensible initial guesses
A = hypot(Ar,Az);
ai = atan2(Az,Ar);
dn = A - B.Geometry.A0;
Q0 = hertz_contactlaw(B.Contact.K,B.Contact.n,dn,B.Contact.tol);
wi = Q0.*cos(ai)/Kri;
wo = Q0.*cos(ai)/Kro;

%make everything 1D
sz = size(Az);
Ar = permute(Ar(:),[2 3 1]);
Az = permute(Az(:),[2 3 1]);
wi = permute(wi(:),[2 3 1]);
wo = permute(wo(:),[2 3 1]);

x = [wi; wo];
dQ = x + Inf;
iter = 0;

iFlex = [];
if B.Options.bRaceCompliancei
    iFlex(end+1) = 1;
end
if B.Options.bRaceComplianceo
    iFlex(end+1) = 2;
end

Qri = Ar + NaN; Kri = Ar + NaN;
Qro = Ar + NaN; Kro = Ar + NaN;

Ar0 = Ar;

while any(abs(dQ(:))>1E-8)
    wi = x(1,:,:);
    wo = x(2,:,:);
    
    Ar = Ar0 - wi - wo;
    A = hypot(Ar,Az);
    ai = atan2(Az,Ar);
    ao = ai;

    dn = A - B.Geometry.A0;
    
    [Qi,Ki] = hertz_contactlaw(B.Contact.K,B.Contact.n,dn,B.Contact.tol);
    Qo = Qi; Ko = Ki;

    if B.Options.bRaceCompliancei
        [Qri,Kri] = race_compliance_loads(B.Race.Inner,-wi); Qri = -Qri;
    end
    if B.Options.bRaceComplianceo
        [Qro,Kro] = race_compliance_loads(B.Race.Outer,wo);
    end

    Kmat = [Ki+Kri Ki;
             Ko   Ko+Kro];

    dQ = [Qi.*cos(ai)-Qri;
          Qo.*cos(ao)-Qro];  

    x(iFlex,:,:) = x(iFlex,:,:) + mtimesx(minvx(Kmat(iFlex,iFlex,:)),dQ(iFlex,:,:));
    iter = iter + 1;
    if iter > 100
        break
    end
end

wi = reshape(x(1,:,:),sz);
wo = reshape(x(2,:,:),sz);