function [F,V,S] = linear_model(B, States) 
xBearing = [States.qi; States.qo];
dxBearing = [States.qidot; States.qodot];
ddxBearing = [States.qiddot; States.qoddot];

fb = B.Kb*xBearing + B.Cb*dxBearing + B.Mb*ddxBearing;

NPts = size(xBearing,2);
F.F = fb;
F.FInt     = zeros(0,NPts);
F.xInt     = zeros(0,NPts);
F.xdotInt  = zeros(0,NPts);
F.xddotInt = zeros(0,NPts);

V = struct();

if nargout > 2
    S.K = repmat(B.Kb,1,1,NPts);
    S.C = repmat(B.Cb,1,1,NPts);
    S.M = repmat(B.Mb,1,1,NPts);
    
    S.Kqq = S.K;
    S.Cqq = S.C;
    S.Mqq = S.M;
    
    S.Kqx = zeros(8,0,NPts);
    S.Kxq = zeros(0,8,NPts);
    S.Kxx = zeros(0,0,NPts);

    S.Cqx = zeros(8,0,NPts);
    S.Cxq = zeros(0,8,NPts);
    S.Cxx = zeros(0,0,NPts);
    
    S.Mqx = zeros(8,0,NPts);
    S.Mxq = zeros(0,8,NPts);
    S.Mxx = zeros(0,0,NPts);
end
