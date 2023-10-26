function [F,V,S] = linear_model(Params, States) 
xBearing = States.qi - States.qo;
dxBearing = States.qidot- States.qodot;
ddxBearing = States.qiddot- States.qoddot;

fb = Params.F0 + Params.K*xBearing + Params.C*dxBearing + Params.M*ddxBearing;

NPts = size(xBearing,2);
F.F = [fb; -fb];

V = struct();

if nargout > 2
    K = [Params.K -Params.K; -Params.K Params.K];
    C = [Params.C -Params.C; -Params.C Params.C];
    M = [Params.M -Params.M; -Params.M Params.M];

    S.K = repmat(K,1,1,NPts);
    S.C = repmat(C,1,1,NPts);
    S.M = repmat(M,1,1,NPts);
end