function [k,m,wFit] = rotor_tune(P,Kt,k0,Mt,m0,iModes,wTarget)

x0 = [k0; m0];

K0 = P.Model.K;
M0 = P.Model.M;

[V,omega] = get_frequencies(K0,M0);
if isempty(iModes)
    ii = get_matching_modes(omega,wTarget);
else
    ii = iModes;
end
omega = omega(ii);
omega2 = sqrt(diag((V(:,ii)'*K0*V(:,ii)))./diag((V(:,ii)'*M0*V(:,ii))));

fun = @(x,auxdata)obj_fun(x,auxdata);
grad = @(x,auxdata)grad_fun(x,auxdata);

%setup problem for ipopt
options.cl = [];
options.cu = [];
options.lb = 0*x0+1E-3;
options.ub = Inf*x0;
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.print_level = 5;
options.ipopt.max_iter = 30;
options.auxdata = {x0,K0,M0,Kt,Mt,V(:,ii),iModes,wTarget};

funcs.objective = fun;
funcs.gradient = grad;
funcs.constraints = [];
funcs.jacobian = [];
funcs.iterfunc = @iteration;

funcs.jacobianstructure = [];

xi = 0*x0 + 1;
% xi = x0;
err = obj_fun(xi,options.auxdata);
for i = 1%:3
    [x, info] = ipopt_auxdata(xi,funcs,options);
    xi = x;
end
x = x0.*x;

[K,M] = get_matrices(x,x0,K0,M0,Kt,Mt);
[~,omega] = get_frequencies(K,M);
if isempty(iModes)
    ii = get_matching_modes(omega,wTarget);
else
    ii = iModes;
end
wFit = omega(ii);

k = x(1:length(Kt));
m = x(length(Kt)+1:end);

function err = obj_fun(x,auxdata)
[x0,K0,M0,Kt,Mt,V,iModes,omegaTarget] = deal(auxdata{:});
x = x0.*x;
[K,M] = get_matrices(x,x0,K0,M0,Kt,Mt);
[~,omega] = get_frequencies(K,M);
if isempty(iModes)
    iModes = get_matching_modes(omega,omegaTarget);
end
omega = omega(iModes);
% omega = sqrt(diag((V'*K*V))./diag((V'*M*V)));
err = sum(((omega - omegaTarget)).^2);
if iscomplex(err)
    keyboard
end

function iModes = get_matching_modes(omega,omegaTarget)
w_fit = omega;
for i = 1:length(omegaTarget)
    w_err = w_fit-omegaTarget(i);
    [~,iModes(i)] = min(w_err.^2,[],'omitnan'); 
    w_fit(iModes(i)) = NaN;
end

function [K,M] = get_matrices(x,x0,K0,M0,Kt,Mt)
dx = x - x0;
dk = dx(1:length(Kt));
dm = dx(length(Kt)+1:end);

K = K0;
for i = 1:length(Kt)
    K = K + dk(i)*Kt{i};
end

M = M0;
for i = 1:length(Mt)
    M = M + dm(i)*Mt{i};
end

function [V,D] = get_frequencies(K,M)
[V,D] = eig(K,M,'vector');
[D,ii] = sort(sqrt(D));
D = real(D);
V = V(:,ii);

function G = grad_fun(x,auxdata)
G = jacobian(@obj_fun,x,auxdata);
% [x0,K0,M0,Kt,Mt,iModes,omegaTarget] = deal(auxdata{:});
% x = x.*x0;