function [F,V,S] = radial_model(B, States) 
xBearing = States.qi - States.qo;
dxBearing = States.qidot - States.qodot;
ddxBearing = States.qiddot - States.qoddot;


x = xBearing([1 3],:);
xdot = dxBearing([1 3],:);
xddot = ddxBearing([1 3],:);

r = sqrt(sum(x.^2));
n = x ./ (repmat(r,2,1) + eps);
t = [n(2,:); -n(1,:)];

rdot = xdot .* n;
rddot = xddot .* n;

fel = getforces(B,r);
fdamp = B.Crad * rdot;
finer = B.Mrad * rddot;

NPts = size(States.qi,2);
fb = zeros(4,NPts);
fb([1 3],:) = (fel + fdamp + finer).*n;
fb([2,4],:) = B.Krot*xBearing([2,4],:) + B.Crot*dxBearing([2,4],:) + B.Mrot*ddxBearing([2,4],:);

F.F = [fb; -fb];
V.Fr = fb;
V.r  = r;
V.rdot = rdot;

if nargout > 2
    [fel,Kr] = getforces(B,r);

    Krad = zeros(2,2,NPts);
    Krad(1,1,:) = Kr;
    Krad(2,2,:) = fel./(r+eps);

    Crad = zeros(2,2,NPts);
    Crad(1,1,:) = B.C;
    
    Mrad = zeros(2,2,NPts);
    Mrad(1,1,:) = B.M;

    R(:,1,:) = permute(n,[1 3 2]);
    R(:,2,:) = permute(t,[1 3 2]);
    Rt = mtransposex(R);

    Kb = zeros(4,4,size(xBearing,2));
    Kb([1 3],[1 3],:) = mtimesx(Rt,mtimesx(Krad,R));
    Kb([2 4],[2 4],:) = B.Krot;

    Cb = zeros(4,4,size(xBearing,2));
    Cb([1 3],[1 3],:) = mtimesx(Rt,mtimesx(Crad,R));
    Cb([2 4],[2 4],:) = B.Crot;

    Mb = zeros(4,4,size(xBearing,2));
    Mb([1 3],[1 3],:) = mtimesx(Rt,mtimesx(Mrad,R));
    Mb([2 4],[2 4],:) = B.Mrot;

    S.K = [Kb -Kb; -Kb Kb];
    S.C = [Cb -Cb; -Cb Cb];
    S.M = [Mb -Mb; -Mb Mb];
end

function [fel,Kr] = getforces(B,r)
fel = B.K * max(r-B.c,0).^B.n;

if nargout > 1
    Kr = B.K * B.n * max(r-B.c,0).^(B.n-1) .* (r > B.c);
end