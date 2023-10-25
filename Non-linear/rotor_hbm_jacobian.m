function Jx = rotor_hbm_jacobian(dfhbm_dxhbm,hbm,problem,w0)

NComp = hbm.harm.NComp;
Nfft = hbm.harm.Nfft;

P = problem.P;
REB = P.Bearing{1}.Params{2};
iRot = P.Model.iRot;
Z = REB.Setup.Z;

NDof = P.Model.NDof;
ijacobx = hbm.nonlin.hbm.ijacobx;

if P.Model.bCompressREB
    K = dfhbm_dxhbm;
    K(1:NDof,NDof+1:end,:) = 0; %remove Kqx
    J = hbm.nonlin.hbm.Jx.*K(ijacobx,ijacobx,:);

    Kqx = 0*dfhbm_dxhbm;
    Kqx(1:NDof,NDof+1:end,:) = dfhbm_dxhbm(1:NDof,NDof+1:end,:);
    Jqxi0 = hbm.nonlin.hbm.Jifft.*Kqx(ijacobx,ijacobx,:);
    Jqx = hbm.nonlin.hbm.Jfft.*shift_balls_fast(Jqxi0,REB,iRot,hbm);
    Jx = J + Jqx;
else
    %shouldn't really use this, but for completeness
    Jx = sum(hbm.nonlin.hbm.Jx.*dfhbm_dxhbm(ijacobx,ijacobx),3);
end
Jx = sum(Jx,3);

function J = shift_balls_fast(J0,REB,iRot,hbm)
J = 0*J0;

if REB.Model.NDof > 0
    Z = REB.Setup.Z;
    for i = 1:Z
        if iRot == 1
           kShift = -(i-1)/Z * hbm.harm.Nfft(1); 
        elseif iRot == 2
            kShift = -(i-1)/Z * hbm.harm.Nfft(2) *  hbm.harm.Nfft(1);
        end

        J = J + circshift(J0,kShift,3);
    end
end

function J = shift_balls(J0,REB,iRot,hbm)
Z = REB.Setup.Z;
Nfft = hbm.harm.Nfft;
J = 0*J0;

J0_3d = reshape(permute(J0,[3 1 2]),Nfft(1),Nfft(2),size(J0,1),size(J0,2));

if REB.Model.NDof > 0
    for i = 1:Z
        kShift = -(i-1)/Z * hbm.harm.Nfft(iRot);
        Ji = circshift(J0_3d,kShift,iRot);
        Ji = ipermute(reshape(Ji,prod(Nfft),size(J,1),size(J,2)),[3 1 2]);
        J = J + Ji;
    end
end