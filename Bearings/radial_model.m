function [F,V,S] = radial_model(B, States) 
xBearing = States.qi - States.qo;
dxBearing = States.qidot - States.qodot;

n = xBearing([1 3],:);
r = sqrt(sum(n.^2));
n = n ./ (repmat(r,2,1) + eps);
rdot = sum(dxBearing([1 3],:) .* n,1);

fr = B.K * max(r-B.c,0).^B.n + B.C * rdot;

fb = zeros(4,size(xBearing,2));
fb([1 3],:) = repmat(fr,[2 1]).*n;

F.F = fb;
F.FInt     = zeros(0,NPts);
F.xInt     = zeros(0,NPts);
F.xdotInt  = zeros(0,NPts);
F.xddotInt = zeros(0,NPts);

V.Fr = fr;
V.r  = r;
V.rdot = rdot;

if nargout > 2
    Kr = B.K * B.n * max(r-B.c,0).^(B.n-1);
    Cr = B.C * (0*r+1);

    n  = permute(n,[1 3 2]);
    nt = permute(n,[2 1 3]);
    N = mtimesx(n,nt);

    Kb = zeros(4,4,size(xBearing,2));
    Kb([1 3],[1 3],:) = repmat(permute(Kr,[1 3 2]),2)*N;

    Cb = zeros(4,4,size(xBearing,2));
    Cb([1 3],[1 3],:) = repmat(permute(Cr,[1 3 2]),2)*N;
    
    S.K = Kb;
    S.C = Cb;
        
    S.Kqq = Kb;
    S.Cqq = Kb;
    
    S.Kqx = zeros(4,0,NPts);
    S.Kxq = zeros(0,4,NPts);
    S.Kxx = zeros(0,0,NPts);

    S.Cqx = zeros(4,0,NPts);
    S.Cxq = zeros(0,4,NPts);
    S.Cxx = zeros(0,0,NPts);
end
