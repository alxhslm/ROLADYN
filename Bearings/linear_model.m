function [F,V,S] = linear_model(B, States) 
xBearing = [States.qi; States.qo];
dxBearing = [States.qidot; States.qodot];
ddxBearing = [States.qiddot; States.qoddot];

fb = B.F0 + B.Kb*xBearing + B.Cb*dxBearing + B.Mb*ddxBearing;

NPts = size(xBearing,2);
F.F = fb;

V = struct();

if nargout > 2
    S.K = repmat(B.Kb,1,1,NPts);
    S.C = repmat(B.Cb,1,1,NPts);
    S.M = repmat(B.Mb,1,1,NPts);
end
