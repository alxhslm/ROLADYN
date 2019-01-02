function Qi = race_contact(Contact,Race,Options,dbi,tol)
if nargin < 3
    tol = 0;
end

% Ki = Contact.K;
% n = Contact.n;
% c = Contact.c;
% e = Contact.e;

% Qi = Ki*sgn_power(dbi,n).*(1+c*abs(dbi).^e);

Qi = ball_newton(dbi,Contact,Race,Options,tol);
% Qi = maxSmooth(Qi,0,tol);

function [Q,iter] = ball_newton(r,Contact,Race,Options,tol)
sz = size(r);
r = permute(r(:),[2 3 1]);
w = [0*r; 0*r];
dQ = w + Inf;
iter = 0;

iFlex = [];
if Options.bRaceCompliancei
    iFlex(end+1) = 1;
end
if Options.bRaceComplianceo
    iFlex(end+1) = 2;
end

Qri = r + NaN; Kri = r + NaN;
Qro = r + NaN; Kro = r + NaN;

while any(abs(dQ(:))>1E-8)
    [Qi,Ki] = hertz_contact(Contact.K,Contact.n,r-w(1,:,:)-w(2,:,:),tol);
    Qo = Qi; Ko = Ki;
    
    if Options.bRaceCompliancei
        [Qri,Kri] = race_compliance(Race.Inner,-w(1,:,:)); Qri = -Qri; %Kri = -Kri;
    end
    if Options.bRaceComplianceo
        [Qro,Kro] = race_compliance(Race.Outer,w(2,:,:));
    end
    
    Kmat = [Ki+Kri Ki;
             Ko   Ko+Kro];

    dQ = [Qi-Qri;
          Qo-Qro];  

    w(iFlex,:,:) = w(iFlex,:,:) + mtimesx(minvx(Kmat(iFlex,iFlex,:)),dQ(iFlex,:,:));
    iter = iter + 1;
    if iter > 100
        Qi = Qi + NaN;
        break
    end
end
Q = reshape(Qi,sz);

function I = minvx(M)
if size(M,1) == 1
    I = 1./M;
elseif size(M,1) == 2
    det = M(1,1,:).*M(2,2,:) - M(1,2,:).*M(2,1,:);
    I = [M(2,2,:) -M(1,2,:);
        -M(2,1,:)  M(1,1,:)]./det;
end