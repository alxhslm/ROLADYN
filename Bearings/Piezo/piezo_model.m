function [Forces,Channels,Stiffness] = piezo_model(Params, States)
%convert bearing displacememnts into the q vector
if size(States.qi,1) == 4
    R = [1     0     0     0];
else
    R = 1;
end

model = Params.fun;
NPts = size(States.qi,2);

States = default_inputs(States,R);
States = default_int(States,NPts);

if nargout < 3
    [Forces,Channels] = feval(model,Params,States);
    Forces.F  = mkron([1;-1],R'*Forces.F);
else
    [Forces,Channels,Stiffness] = stiffnessAndDamping(model,Params,States,R);
end

Forces.xInt     = States.xInt;
Forces.xdotInt  = States.xdotInt;
Forces.xddotInt = States.xddotInt;

Forces.q = R'*(States.qi - States.qo);

function [Forces,Channels,Stiffness] = stiffnessAndDamping(model,Params,States,R)
[Forces,Channels,Stiffness] = feval(model,Params,States);

if ~isempty(Stiffness)
    h = 1E-8;
    q0 = States.qi;
    qdot0 = States.qidot;
    qddot0 = States.qiddot;
    xInt0 = States.xInt;
    xIntdot0 = States.xIntdot;
    xIntddot0 = States.xIntddot;
    F0 = Forces.F;
    FInt0 = Forces.FInt;

    %stiffness
    States.qi = q0 + h;
    Forces = feval(model,Params,States);
    States.qi = q0; 
    Stiffness.Kqq = permute((Forces.F-F0)./h,[1 3 2]);
    Stiffness.Kxq = permute((Forces.FInt-FInt0)./h,[1 3 2]);

    States.xInt = xInt0 + h;
    Forces = feval(model,Params,States);
    States.xInt = xInt0; 
    Stiffness.Kqx = permute((Forces.F-F0)./h,[1 3 2]);
    Stiffness.Kxx = permute((Forces.FInt-FInt0)./h,[1 3 2]);

    %damping 
    States.qidot = qdot0 + h;
    Forces = feval(model,Params,States);
    States.qidot = qdot0; 
    Stiffness.Cqq = permute((Forces.F-F0)./h,[1 3 2]);
    Stiffness.Cxq = permute((Forces.FInt-FInt0)./h,[1 3 2]);

    States.xIntdot = xIntdot0 + h;
    Forces = feval(model,Params,States);
    States.xIntdot = xIntdot0; 
    Stiffness.Cqx = permute((Forces.F-F0)./h,[1 3 2]);
    Stiffness.Cxx = permute((Forces.FInt-FInt0)./h,[1 3 2]);

    %inertia 
    States.qiddot = qddot0 + h;
    Forces = feval(model,Params,States);
    States.qiddot = qddot0; 
    Stiffness.Mqq = permute((Forces.F-F0)./h,[1 3 2]);
    Stiffness.Mxq = permute((Forces.FInt-FInt0)./h,[1 3 2]);

    States.xIntddot = xIntddot0 + h;
    Forces = feval(model,Params,States);
    States.xIntddot = xIntddot0; 
    Stiffness.Mqx = permute((Forces.F-F0)./h,[1 3 2]);
    Stiffness.Mxx = permute((Forces.FInt-FInt0)./h,[1 3 2]);
end

Forces.F  = mkron([1;-1],R'*Forces.F);

%Combine stiffness terms
Stiffness.K   = mkron([1 -1; -1 1],mtimesx(R',mtimesx(Stiffness.K,R)));
Stiffness.Kqq = mkron([1 -1; -1 1],mtimesx(R',mtimesx(Stiffness.Kqq,R)));
Stiffness.Kqx = mkron([1; -1],mtimesx(R',Stiffness.Kqx));
Stiffness.Kxq = mkron([1 -1],mtimesx(Stiffness.Kxq,R));

%Combine damping terms
Stiffness.C   = mkron([1 -1; -1 1],mtimesx(R',mtimesx(Stiffness.C,R)));
Stiffness.Cqq = mkron([1 -1; -1 1],mtimesx(R',mtimesx(Stiffness.Cqq,R)));
Stiffness.Cqx = mkron([1; -1],mtimesx(R',Stiffness.Cqx));
Stiffness.Cxq = mkron([1 -1],mtimesx(Stiffness.Cxq,R));

%Combine inertia terms
Stiffness.M   = mkron([1 -1; -1 1],mtimesx(R',mtimesx(Stiffness.M,R)));
Stiffness.Mqq = mkron([1 -1; -1 1],mtimesx(R',mtimesx(Stiffness.Mqq,R)));
Stiffness.Mqx = mkron([1; -1],mtimesx(R',Stiffness.Mqx));
Stiffness.Mxq = mkron([1 -1],mtimesx(Stiffness.Mxq,R));

function States = default_int(States,NPts)
if ~isfield(States,'bSolve')
    States.bSolve = 1;
end
if ~isfield(States,'xInt')
    States.xInt = zeros(0,NPts);
end
if ~isfield(States,'xIntdot')
    States.xIntdot = 0*States.xInt;
end
if ~isfield(States,'xIntddot')
    States.xIntddot = 0*States.xInt;
end

function States = default_inputs(States,R)
if ~isfield(States,'qo')
    States.qo = 0*States.qi;
end
suffix = {'dot','ddot'};
IO = {'i','o'};
for i = 1:length(suffix)
    for k = 1:2
        States.(['q' IO{k} suffix{i}]) = 0*States.qi;
    end
end

prefix = {'q','q','q'};
suffix = {'' ,'dot','ddot'};
IO = {'i','o'};
for i = 1:length(prefix)
    for k = 1:2
        States.([prefix{i} IO{k} suffix{i}]) = R*States.([prefix{i} IO{k} suffix{i}]);
    end
end