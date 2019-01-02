function [F,V,S] = empty_model(States)
NPts = size(States.qi,2);

F.F = zeros(8,NPts);
F.FInt     = zeros(0,NPts);
F.xInt     = zeros(0,NPts);
F.xdotInt  = zeros(0,NPts);
F.xddotInt = zeros(0,NPts);

V = struct();

if nargout > 2
    S.K = zeros(8,8,NPts);
    S.C = zeros(8,8,NPts);
    
    S.Kqq = S.K;
    S.Cqq = S.C;
    
    S.Kqx = zeros(8,0,NPts);
    S.Kxq = zeros(0,8,NPts);
    S.Kxx = zeros(0,0,NPts);

    S.Cqx = zeros(8,0,NPts);
    S.Cxq = zeros(0,8,NPts);
    S.Cxx = zeros(0,0,NPts);
end
