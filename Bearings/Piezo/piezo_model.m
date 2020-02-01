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

if isempty(Stiffness)
    Stiffness = struct();
    
    h = 1E-10;
    F0 = Forces.F;
    FInt0 = Forces.FInt;
    NDofInt = size(States.xInt,1);
    NPts = size(States.qi,2);
        
    %stiffness
    q0 = States.qi;
    States.qi = q0 + h;
    Forces = feval(model,Params,States);
    Stiffness.Kqq = permute((Forces.F-F0)./h,[1 3 2]);
    Stiffness.Kxq = permute((Forces.FInt-FInt0)./h,[1 3 2]);
    States.qi = q0; 

    xInt0 = States.xInt;
    Stiffness.Kqx = zeros(1,NDofInt,NPts);
    Stiffness.Kxx = zeros(NDofInt,NDofInt,NPts);
    for i = 1:NDofInt
        States.xInt(i,:) = xInt0(i,:) + h;
        Forces = feval(model,Params,States);
        Stiffness.Kqx(:,i,:) = permute((Forces.F-F0)./h,[1 3 2]);
        Stiffness.Kxx(:,i,:) = permute((Forces.FInt-FInt0)./h,[1 3 2]);
        States.xInt = xInt0; 
    end

    %damping 
    qdot0 = States.qidot;
    States.qidot = qdot0 + h;
    Forces = feval(model,Params,States);
    Stiffness.Cqq = permute((Forces.F-F0)./h,[1 3 2]);
    Stiffness.Cxq = permute((Forces.FInt-FInt0)./h,[1 3 2]);
    States.qidot = qdot0; 
    
    xdotInt0 = States.xdotInt;
    Stiffness.Cqx = zeros(1,NDofInt,NPts);
    Stiffness.Cxx = zeros(NDofInt,NDofInt,NPts);
    for i = 1:NDofInt
        States.xdotInt(i,:) = xdotInt0(i,:) + h;
        Forces = feval(model,Params,States);
        Stiffness.Cqx(:,i,:) = permute((Forces.F-F0)./h,[1 3 2]);
        Stiffness.Cxx(:,i,:) = permute((Forces.FInt-FInt0)./h,[1 3 2]);
        States.xdotInt = xdotInt0; 
    end
    
    %input
    u0 = States.u;
    States.u = u0 + h;
    Forces = feval(model,Params,States);
    Stiffness.Kqu = permute((Forces.F-F0)./h,[1 3 2]);
    Stiffness.Kxu = permute((Forces.FInt-FInt0)./h,[1 3 2]);
    States.u = u0; 

    udot0 = States.udot;
    States.udot = udot0 + h;
    Forces = feval(model,Params,States);
    Stiffness.Cqu = permute((Forces.F-F0)./h,[1 3 2]);
    Stiffness.Cxu = permute((Forces.FInt-FInt0)./h,[1 3 2]);
    States.udot = udot0; 
end

Forces.F  = mkron([1;-1],R'*Forces.F);


if Params.NDofTot > 0
    A = -mtimesx(minvx(Stiffness.Kxx + 1i*Stiffness.Cxx),Stiffness.Kxq + 1i*Stiffness.Cxq);
    B = -mtimesx(minvx(Stiffness.Kxx + 1i*Stiffness.Cxx),Stiffness.Kxu + 1i*Stiffness.Cxu);
    
    Stiffness.K = Stiffness.Kqq + mtimesx(Stiffness.Kqx,A);
    Stiffness.C = Stiffness.Cqq + mtimesx(Stiffness.Cqx,A);

    Stiffness.Ku = Stiffness.Kqu + mtimesx(Stiffness.Kqx,B);
    Stiffness.Cu = Stiffness.Cqu + mtimesx(Stiffness.Cqx,B);
else
    Stiffness.K = Stiffness.Kqq;
    Stiffness.C = Stiffness.Cqq;
    
    Stiffness.Ku = Stiffness.Kqu;
    Stiffness.Cu = Stiffness.Cqu;
        
    Npts = size(States.qi,2);
    Stiffness.Kqx = zeros(1,0,Npts);
    Stiffness.Kxq = zeros(0,1,Npts);
    Stiffness.Kxx = zeros(0,0,Npts);
    
    Stiffness.Cqx = zeros(1,0,Npts);
    Stiffness.Cxq = zeros(0,1,Npts);
    Stiffness.Cxx = zeros(0,0,Npts);
    
    Stiffness.Kxu = zeros(0,1,Npts);
    Stiffness.Cxu = zeros(0,1,Npts);
end


%Combine stiffness terms
Stiffness.Kqq = mkron([1 -1; -1 1],mtimesx(R',mtimesx(Stiffness.Kqq,R)));
Stiffness.Kqx = mkron([1; -1],mtimesx(R',Stiffness.Kqx));
Stiffness.Kxq = mkron([1 -1],mtimesx(Stiffness.Kxq,R));
Stiffness.K   = mkron([1 -1; -1 1],mtimesx(R',mtimesx(Stiffness.K,R)));
Stiffness.Ku  = mkron([1; -1],mtimesx(R',Stiffness.Ku));
Stiffness.Kqu = mkron([1; -1],mtimesx(R',Stiffness.Kqu));

%Combine damping terms
Stiffness.Cqq = mkron([1 -1; -1 1],mtimesx(R',mtimesx(Stiffness.Cqq,R)));
Stiffness.Cqx = mkron([1; -1],mtimesx(R',Stiffness.Cqx));
Stiffness.Cxq = mkron([1 -1],mtimesx(Stiffness.Cxq,R));
Stiffness.C   = mkron([1 -1; -1 1],mtimesx(R',mtimesx(Stiffness.C,R)));
Stiffness.Cu  = mkron([1; -1],mtimesx(R',Stiffness.Cu));
Stiffness.Cqu = mkron([1; -1],mtimesx(R',Stiffness.Cqu));

function States = default_int(States,NPts)
if ~isfield(States,'bSolve')
    States.bSolve = 1;
end
if ~isfield(States,'xInt')
    States.xInt = zeros(0,NPts);
end
if ~isfield(States,'xdotInt')
    States.xdotInt = 0*States.xInt;
end
if ~isfield(States,'xddotInt')
    States.xddotInt = 0*States.xInt;
end

function States = default_inputs(States,R)
if ~isfield(States,'qo')
    States.qo = 0*States.qi;
end

suffix = {'dot','ddot'};
IO = {'i','o'};
for i = 1:length(suffix)
    for k = 1:2
        if ~isfield(States,['q' IO{k} suffix{i}])
            States.(['q' IO{k} suffix{i}]) = 0*States.qi;
        end
    end
end

suffix = {'' ,'dot','ddot'};
IO = {'i','o'};
for i = 1:length(suffix)
    for k = 1:2
        States.(['q' IO{k} suffix{i}]) = R*States.(['q' IO{k} suffix{i}]);
    end
end


if ~isfield(States,'u') || isempty(States.u)
    States.u = repmat(0*States.qi(1,:),2,1);
end

suffix = {'dot','ddot'};
for i = 1:length(suffix)
    if ~isfield(States,['u' suffix{i}])
        States.(['u' suffix{i}]) = 0*States.u;
    end
end