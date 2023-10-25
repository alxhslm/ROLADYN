function [F,V,S] = REB_preload(Params,States,q0)

%decide if we need to show the waitbar
if ~isfield(States,'bWaitbar')
    States.bWaitbar = 0;
end

%solve for the ball equilibrium
NPts = size(States.qi,2);

States = default_inputs(States);
States = default_speeds(States,NPts);

%find free and fixed dof
iConstr = isinf(diag(Params.Setup.KbParallel));
iFree = isnan(States.qi(:,1));
iFixed = isnan(States.F(:,1));
if any(iConstr & iFree)
    error(error('Problem specified incorrectly'));
end
if any(iFixed & iFree)
    error('Problem specified incorrectly')
end
States.qi(iConstr,:) = 0;
N = sum(iFree);

Jstr = ones(N);
x0 = zeros(Params.Model.NDof,1);
if nargin < 3
    q0 = zeros(4,1);
end
X0 = [q0(iFree); x0];

options.print_level = 0;
options.maxit = 400;
    
if States.bWaitbar
    h = waitbar(0,'Sweep');
end
Npts = size(States.qi,2);
for j = 1:Npts
    States_j = getjthpoint(States,j);
    options.jacobian = @preload_jac;
    options.jacobianstructure = Jstr;
    info.status = 2;
    count = 0;
    bConverged = 0;
    while count < 10
        [X,info] = fipopt([],X0,@preload_constr,options,Params,States_j,N,iFree);
        F = preload_constr(X,Params,States_j,N,iFree);
        X0 = X;
        count = count + 1;
        if any(info.status==[0 1])
           bConverged = 1;
        end
    end
    if bConverged
        Xsol(:,j) = X;
    else
        Xsol(:,j) = NaN;
    end
    if States.bWaitbar
        waitbar(j/Npts,h);
    end
end
if States.bWaitbar
    close(h);
end

States.qi(iFree,:) = Xsol(1:N,:);
if nargout > 2
    [F,V,S] = REB_model(Params, States);
else
    [F,V] = REB_model(Params, States);
end
F.q = States.qi;

function States = getjthpoint(States,j)
prefix = {'q','q','q','O','A'};
suffix = {'' ,'dot','ddot','',''};
IO = {'i','o'};
for i = 1:length(prefix)
    for k = 1:2
        States.([prefix{i} IO{k} suffix{i}]) = States.([prefix{i} IO{k} suffix{i}])(:,j);
    end
end

States.F = States.F(:,j);

function States = default_inputs(States)
if ~isfield(States,'qo')
    States.qo = zeros(4,size(States.qi,2));
end
suffix = {'dot','ddot'};
IO = {'i','o'};
for i = 1:length(suffix)
    for k = 1:2
        States.(['q' IO{k} suffix{i}]) = zeros(4,size(States.qi,2));
    end
end

prefix = {'q','q','q'};
suffix = {'' ,'dot','ddot'};
IO = {'i','o'};
for i = 1:length(prefix)
    for k = 1:2
        States.([prefix{i} IO{k} suffix{i}]) = States.([prefix{i} IO{k} suffix{i}]);
    end
end

function States = default_speeds(States,NPts)
fields = {'Oi','Oo','Ai','Ao'};
for i = 1:length(fields)
    if ~isfield(States,fields{i})
        States.(fields{i}) = 0;
    end
end

if length(States.Oi) == 1
    States.Oi = States.Oi + zeros(1,NPts);
    States.Oo = States.Oo + zeros(1,NPts);
    States.Ai = States.Ai + zeros(1,NPts);
    States.Ao = States.Ao + zeros(1,NPts);
end

function F = preload_constr(X,Params,States,N,iFree)
States.qi(iFree) = X(1:N);
Forces = REB_model(Params,States);
FFree = Forces.F(1:4,:) - States.F;
F = FFree(iFree);

function J = preload_jac(X,Params,States,N,iFree)
States.qi(iFree) = X(1:N);
[F,~,S] = REB_model(Params,States);
J = S.Kqq(iFree,iFree);