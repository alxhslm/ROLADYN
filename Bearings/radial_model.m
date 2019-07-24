function [F,V,S] = radial_model(B, States) 
xBearing = States.qi - States.qo;
dxBearing = States.qidot - States.qodot;

x = xBearing([1 3],:);
r = sqrt(sum(x.^2));
n = x ./ (repmat(r,2,1) + eps);
t = [-n(2,:); n(1,:)];

rdot = sum(dxBearing([1 3],:) .* n,1);

tol = 1E-16;
if B.c > 0 && tol > 0
    fr = B.K * maxSmooth(r-B.c,0,tol).^B.n + B.C * rdot;
else
    fr = B.K * max(r-B.c,0).^B.n + B.C * rdot;
end
    
fb = zeros(4,size(xBearing,2));
fb([1 3],:) = fr.*n;

F.Fi =  fb;
F.Fo = -fb;
F.F = [F.Fi; F.Fo];
NPts = size(States.qi,2);
F.FInt     = zeros(0,NPts);
F.xInt     = zeros(0,NPts);
F.xdotInt  = zeros(0,NPts);
F.xddotInt = zeros(0,NPts);

V.Fr = fr;
V.r  = r;
V.rdot = rdot;

if nargout > 2
    if B.c > 0 && tol > 0
        [rContact,drContact] = maxSmooth(r-B.c,0,tol);
        Kr = B.K * (B.n * rContact.^(B.n-1) + rContact.^B.n.*drContact);
    else
        Kr = B.K * B.n * max(r-B.c,0).^(B.n-1);
    end
    Cr = B.C * (0*r+1);
    
    z = zeros(1,1,NPts);
    Krot = [permute(Kr,[1 3 2]) z;
                  z      permute(fr./(r+eps),[1 3 2])];
        
    Crot = [permute(Cr,[1 3 2]) z;
            z  z];

    R(:,1,:) = permute(n,[1 3 2]);
    R(:,2,:) = permute(t,[1 3 2]);
    Rt = mtransposex(R);

    Kb = zeros(4,4,size(xBearing,2));
    Kb([1 3],[1 3],:) = mtimesx(Rt,mtimesx(Krot,R));
    
    Cb = zeros(4,4,size(xBearing,2));
    Cb([1 3],[1 3],:) = mtimesx(Rt,mtimesx(Crot,R));
    
    K = [Kb -Kb; -Kb Kb];
    C = [Cb -Cb; -Cb Cb];
        
    S.K = K;
    S.C = C;

    S.Kqq = K;
    S.Cqq = C; 
    
    S.Kqx = zeros(8,0,NPts);
    S.Kxq = zeros(0,8,NPts);
    S.Kxx = zeros(0,0,NPts);

    S.Cqx = zeros(8,0,NPts);
    S.Cxq = zeros(0,8,NPts);
    S.Cxx = zeros(0,0,NPts);
end
