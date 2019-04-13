function [Qi,Qo,wi,wo,K] = race_compliance_contactlaw(Contact,Race,Options,dn)

    [wi,wo] = ball_newton(Contact,Race,Options,dn);
    
    [Qi,Kh] = hertz_contactlaw(Contact.K,Contact.n,dn,Contact.tol);
    Qo = Qi;
    [~,Kri] = race_compliance_loads(Race.Inner,-wi);
    [~,Kro] = race_compliance_loads(Race.Outer,wo);
    
    K = 1./(1./Kh + 1./Kri + 1./Kro);
    
    function [wi,wo,iter] = ball_newton(Contact,Race,Options,dn)
    sz = size(dn);
    dn = permute(dn(:),[2 3 1]);
    w = [0*dn; 0*dn];
    dQ = w + Inf;
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
        [Qi,Ki] = hertz_contactlaw(Contact.K,Contact.n,dn-w(1,:,:)-w(2,:,:),Contact.tol);
        Qo = Qi; Ko = Ki;
        
        if Options.bRaceCompliancei
            [Qri,Kri] = race_compliance_loads(Race.Inner,-w(1,:,:)); Qri = -Qri; %Kri = -Kri;
        end
        if Options.bRaceComplianceo
            [Qro,Kro] = race_compliance_loads(Race.Outer,w(2,:,:));
        end
        
        Kmat = [Ki+Kri Ki;
                 Ko   Ko+Kro];
    
        dQ = [Qi-Qri;
              Qo-Qro];  
    
        w(iFlex,:,:) = w(iFlex,:,:) + mtimesx(minvx(Kmat(iFlex,iFlex,:)),dQ(iFlex,:,:));
        iter = iter + 1;
        if iter > 100
            break
        end
    end
    
    wi = reshape(w(1,:,:),sz);
    wo = reshape(w(2,:,:),sz);