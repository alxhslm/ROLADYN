function [F,V,S] = empty_model(Params, States)
NPts = size(States.qi,2);

F.F = zeros(8,NPts);

V = struct();

if nargout > 2
    S.K = zeros(8,8,NPts);
    S.C = zeros(8,8,NPts);
    S.M = zeros(8,8,NPts);
end
