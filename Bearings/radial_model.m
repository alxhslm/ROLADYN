function [F,V,S] = radial_model(B, States) 
xBearing = States.qi - States.qo;
dxBearing = States.qidot - States.qodot;

x = xBearing([1 3],:);
xdot = dxBearing([1 3],:);

r = sqrt(sum(x.^2));
n = x ./ (repmat(r,2,1) + eps);
t = [n(2,:); -n(1,:)];

rdot = xdot .* n;

fel = getforces(B,r);
fdamp = B.C * rdot;

NPts = size(States.qi,2);
fb = zeros(4,NPts);
fb([1 3],:) = (fel + fdamp).*n;

F.Fi =  fb;
F.Fo = -fb;
F.F = [F.Fi; F.Fo];
F.FInt     = zeros(0,NPts);
F.xInt     = zeros(0,NPts);
F.xdotInt  = zeros(0,NPts);
F.xddotInt = zeros(0,NPts);

V.Fr = fb;
V.r  = r;
V.rdot = rdot;

if nargout > 2
    [fel,Kr] = getforces(B,r);

    Krot = zeros(2,2,NPts);
    Krot(1,1,:) = Kr;
    Krot(2,2,:) = fel./(r+eps);

    Crot = zeros(2,2,NPts);
    Crot(1,1,:) = B.C;

    R(:,1,:) = permute(n,[1 3 2]);
    R(:,2,:) = permute(t,[1 3 2]);
    Rt = mtransposex(R);

    Kb = zeros(4,4,size(xBearing,2));
    Kb([1 3],[1 3],:) = mtimesx(Rt,mtimesx(Krot,R));

    Cb = zeros(4,4,size(xBearing,2));
    Cb([1 3],[1 3],:) = mtimesx(Rt,mtimesx(Crot,R));

    % h = 1E-9;
    % Kb = zeros(4,4,NPts);
    % for i = [1 3]
    %     States.qi(i,:) = States.qi(i,:) + h;
    %     Forces = radial_model(B, States); 
    %     Kb(:,i,:) = permute(Forces.Fi - fb,[1 3 2])/h;
    %     States.qi(i,:) = States.qi(i,:) - h;
    % end
    
    % Cb = zeros(4,4,NPts);
    % for i = [1 3]
    %     States.qidot(i,:) = States.qidot(i,:) + h;
    %     Forces = radial_model(B, States); 
    %     Cb(:,i,:) = permute(Forces.Fi - fb,[1 3 2])/h;
    %     States.qidot(i,:) = States.qidot(i,:) - h;
    % end
    
    
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


function [fel,Kr] = getforces(B,r)
fel = B.K * max(r-B.c,0).^B.n;

if nargout > 1
    Kr = B.K * B.n * max(r-B.c,0).^(B.n-1) .* (r > B.c);
end