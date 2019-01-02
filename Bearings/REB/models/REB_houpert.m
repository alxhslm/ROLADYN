function [F,V,S] = REB_houpert(B,States)

q    = States.qi    - States.qo;
Oi = States.Oi;
Oo = States.Oo;

wons = (B.psi*0+1);
q0 = (q(1,:)*0 + 1);

Dx = q(1,:) + q(5,:)*B.z - q(5,:)*B.Ri*tan(B.alpha0);
Dy = q(2,:) - q(4,:)*B.z + q(4,:)*B.Ri*tan(B.alpha0);
Dr = sqrt(Dx.^2 + Dy.^2) - B.cr;
psiR = atan2(Dy,Dx);
A = (q(3,:) - B.cz)*tan(B.alpha0) ./ (Dr+eps);

if B.bCentrifugal 
    %select the ball speed ratio depending on the operating condition
    Ocage = (1 - B.dm/B.D * cos(B.alpha0))./2 * Oi + (1 + B.dm/B.D * cos(B.alpha0))./2 * Oo;

    %centrifugal loads
    Fcent = 0.5*B.mBall*B.dm*Ocage.^2;
    
    gamma = Fcent/cos(B.alpha0)./(B.Inner.K*sgn_power(Dr*cos(B.alpha0),B.n));
    lambda = (B.Outer.K / B.Inner.K)^(1/B.n);
    
    [La,Lr] = integrals_centrifugal(gamma,lambda,A,B.n);
    
    Fr =  B.Inner.K * ((Dr*cos(B.alpha0)).^B.n) * cos(B.alpha0) * B.Setup.Z * Lr;
    Fz =  B.Inner.K * ((Dr*cos(B.alpha0)).^B.n) * sin(B.alpha0) * B.Setup.Z * La;
    
else
    if abs(Dr)  < 1E-8
        Fz = B.K * max(q(3,:) - B.cz,0).^B.n * sin(B.alpha0) * B.Setup.Z;
        Fr = 0*Fz;
    else
        [La,Lr] = integrals(A,B.n);

        Fr =  B.K * ((Dr*cos(B.alpha0)).^B.n) * cos(B.alpha0) * B.Setup.Z * Lr;
        Fz =  B.K * ((Dr*cos(B.alpha0)).^B.n) * sin(B.alpha0) * B.Setup.Z * La;
    end

end

%just divide the load equally between the balls
Q = sqrt(Fr.^2 + Fz.^2)/B.Setup.Z;

Fx = Fr .* cos(psiR);
Fy = Fr .* sin(psiR);

W = [Fx;
     Fy;
     Fz;
     B.Ri*tan(B.alpha0)*Fy;
    -B.Ri*tan(B.alpha0)*Fx;
     0*q0];

%conact angles are constant by definition
alpha_i = B.alpha*q0;
alpha_o = B.alpha*q0;

Qi = Q;
Qo = Q;

%forces
F.F = W - B.cHydro;
F.FInt = zeros(0,size(W,2));

%contact loads
V.Qi = Qi; V.Qo = Qo; 
V.Fi = Fi; V.Fo = Fo;

%geometry
V.alpha_i = alpha_i; V.alpha_o = alpha_o;

%stiffnesses
S = struct([]);

function [La,Lr] = integrals_centrifugal(gamma,lambda,A,n)

p0a = [0.2716
    2.2174
   -3.8974
   -0.2800
   -0.5885
    1.2328
    0.0303
];

p0r = [0.2581
    1.5358
   -9.5203
    0.0194
   -0.9801
   -0.2838
    0.0161
    ];

La = poly(gamma,lambda,A,p0a,n);
Lr = poly(gamma,lambda,A,p0r,n);

function y = poly(gamma,lambda,A,p,n)
gamma = gamma + eps;
f = exp(-p(7)*gamma);
y = maxSmooth(p(1) * f.* (1 + A).^p(2) .* (1 - 1./(1+lambda)).^n + (1-f).*(p(3).*gamma.^p(4) .* lambda.^p(5)).*(1 + A).^p(6),0,1E-6);


function [La,Lr] = integrals(A,n)
if A <= 1
    La = 0.282*(1+A).^2.032;
    Lr = 0.218*(1+A).^1.90;
else
    La = 0.431*(1+A).^1.776;
    Lr = 0.560*(1+A).^0.611;
end
return;

N = 100;
psi = linspace(-pi, pi,N+1)';

wons = 0*psi + 1;
A0 = A*0 + 1;
integrand = max(cos(psi)*A0 + wons*A,0).^n;
La = 1/(2*pi) * trapz(integrand,1) * 2*pi/N;
Lr = 1/(2*pi) * trapz(integrand.*(cos(psi)*A0),1)*2*pi/N;