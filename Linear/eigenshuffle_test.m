function [V,d,W,p] = eigenshuffle_test

fun= @myfun;
p = [0 4];

p0 = p(1);
pEnd = p(end);

A = feval(fun,p0);
[V,d,W] = eig(A,'vector');
W = conj(W);
[V,W] = mkorthog(A,V,W);

% V = V(:,1);
% W = W(:,1);
% d = d(1);

NDof = size(A,1);
NMode = size(V,2);
x0 = packeig(V,d,W);
[V2,d2,W2] = unpackeig(x0,NDof,NMode);

Xscale = 0*x0 + norm(x0);
Xscale(end+1) = 1;%mean([p0,pEnd]);

pMax = max(pEnd,p0);
pMin = min(pEnd,p0);

fun_constr = @(x)eigen_cont(x,fun,NMode,Xscale);
fun_jacobian = @(x)eigen_jac(x,fun,NMode,Xscale);

p = p0;
x = x0;
Xprev = [x0; p0]./Xscale;
F = fun_constr(Xprev);
J = fun_jacobian(Xprev);

% X = null(J);
% [v,d] = unpackeig(X(1:end-1,1).*Xscale(1:end-1),NDof,NMode);
% [v2,d2] = unpackeig(X(1:end-1,2).*Xscale(1:end-1),NDof,NMode);
% % (A*v)./v
% % 
% F2 = fun_constr(Xprev+1E-6*X(:,1)./Xscale);

step = 0.0001;
max_step = 12;
min_step = 0.00001;
count = 0;
options.xtol = 1E-8;
options.ftol = 1E-8;
options.C = 1.1;
options.c = 0.5;
options.maxit = 100;

hWaitbar = waitbar(0, 'Frequency Range');
J = feval(fun_jacobian,Xprev);
tangent_prev = mynull(J);
tangent_prev = tangent_prev * sign(tangent_prev(end)) * sign(pEnd - p0);

fCamp = figure;
subplot(211)
hR = plot(p0,real(d),'o-');
xlim([pMin pMax]);
ylabel('Real')
subplot(212)
hI = plot(p0,imag(d),'o-');
xlim([pMin pMax]);
xlabel('Speed (rad/s)')
ylabel('Imag')
drawnow

% p = linspace(pMin,pMax,100);
% for i = 1:length(p)
%     p0 = p(i);
%     A = feval(fun,p0);
%     [V,d] = eig(A,'vector');
%     x0 = packeig(V,d);
%     X = [x0; p0]./Xscale;
%     F = fun_constr(X);
% end


while p(end) <= pMax && p(end) >= pMin
    X = Xprev + step*tangent_prev;
    iter = 0;
    Xlast = X + Inf;
    tangent = tangent_prev;
    while (norm(X - Xlast) > options.xtol || norm(F) > options.ftol)
        Xlast = X;
        J = feval(fun_jacobian,X);
        F = feval(fun_constr,X);
        if 1
            X = Xlast - pinv(J)*F;
            tangent = mynull(J);
        else
            B = [J; tangent'];
            R = [J*tangent; 0];
            Q = [F; 0];
            W = tangent - B\R;
            tangent = W/norm(W);
            X = Xlast - B\Q;
        end
        
        iter = iter + 1;
        if iter > options.maxit
             break;
        end
    end
    
    pCurr = X(end).*Xscale(end);
    
    if iter < options.maxit && (pCurr - p(end))*(pEnd - p0) > 0
        step = norm(X - Xprev);
        if step > max_step
            %too large step
            step = step * max_step / step;
        elseif step < min_step
            step = step * min_step / step;
        else
            %just right
            count = 0;
            Xprev = X;

            if iter < 5
                step = min(step * options.C,step * max_step / step);
            elseif iter > 10
                step = step / options.C;
            end
            
            tangent_prev = sign(tangent(end))*tangent*sign(pEnd - p0);
            
            %store the data
            p(end+1) = pCurr;
            x(:,end+1) = X(1:end-1).*Xscale(1:end-1);
            [V(:,:,end+1),d(:,end+1)] = unpackeig(x(:,end),NDof,NMode);
            
            for i = 1:NMode
                set(hR(i),'xdata',p,'ydata',real(d(i,:)));
                set(hI(i),'xdata',p,'ydata',imag(d(i,:)));
            end
            
        end
    else
        %failed
        step = step*options.c;
        count = count + 1;
%         if count > 4
%             break
%         end
    end
    
    %update our progress
    waitbar((pCurr-p0)/(pEnd - p0),hWaitbar,sprintf('Step Size: %f',step));
    drawnow
end
close(hWaitbar)
% close(fCamp)

function [V,W] = mkorthog(B,V,W)
 %ensure we normalise the modes correctly
scale = sqrt(diag(W.' * B * V)).';    
scale(scale == 0) = 1;

%now scale the eigenvectors
p = size(V,1);
V = V./ repmat(scale,p,1);
W = W./ repmat(scale,p,1);

function [A,dAdp] = myfun(p)
M = diag([1 1]);
K = [2 -1; -1 1];
Kp = [1 -1; -1 1];
K = K + Kp*p;

C = eye(2)*1E-3;
I = eye(2);
Z = zeros(2);

A = [-M\C -M\K;
      I Z];

dAdp = [-Z -M\Kp;
        Z Z];
    
% A = [1+p 0;
%      0   2];
% dAdp = [1 0;
%         0 0];

function t = mynull(A)
[U,S,V] = svd(A,'econ');
t = V(:,end);

function x = packeig(V,d,W)
d = [real(d); imag(d)];
V = [real(V); imag(V)];
W = [real(W); imag(W)];
x = [d(:);V(:);W(:)];

function [V,d,W] = unpackeig(x,NDof,NMode)
d = x(1:NMode) + 1i*x(NMode+(1:NMode));
% V = reshape(x((N+1):end),N,N);
V = reshape(x((2*NMode) + (1:2*NMode*NDof)),2*NDof,NMode);
V = V(1:NDof,:) + 1i*V((NDof+1):end,:);

W = reshape(x((2*NMode+2*NMode*NDof) + (1:2*NMode*NDof)),2*NDof,NMode);
W = W(1:NDof,:) + 1i*W((NDof+1):end,:);

function f = eigen_cont(X,fun,NModes,Xscale)
x = X(1:end-1).*Xscale(1:end-1);
p = X(end).*Xscale(end);

A = feval(fun,p);
NDof = size(A,1);

[V,d,W] = unpackeig(x,NDof,NModes);

D = diag(d);
residV = (A*V - V*D);
residW = (W.'*A - D*W.');
residD =  diag(W.'*A*V) - 1;

f = packeig(residV,residD,residW);

function J = complex2comp(J,dim)
if length(dim)>1
    N = dim(2);
    M = size(J,2)/N;
    A = zeros(size(J,1),M*N*2);
    for j = 1:M
        A(:,(j-1)*2*N+(1:N)) = J(:,(j-1)*N+(1:N));
        A(:,(j-1)*2*N+N+(1:N)) = -J(:,(j-1)*N+(1:N))/1i;
    end
    J = A;
end

N = dim(1);
M = size(J,1)/N;
A = zeros(M*N*2,size(J,2));
for i = 1:M
    A((i-1)*2*N+(1:N),:) = real(J((i-1)*N+(1:N),:));
    A((i-1)*2*N+N+(1:N),:) = imag(J((i-1)*N+(1:N),:));
end
J = A;

function J = eigen_jac(X,fun,NModes,Xscale)
J = jacobian(@eigen_cont,X,fun,NModes,Xscale);
return

x = X(1:end-1).*Xscale(1:end-1);
p = X(end).*Xscale(end);

[A,dAdp] = feval(fun,p);
NDof = size(A,1);

[V,d,W] = unpackeig(x,NDof,NModes);

D = diag(d);
residD = diag(W.'*A*V) - 1;
residV = (A*V - V*D);
residW = (W.'*A - D*W.');

%jacob of W'*A*V - 1 = 0
Jdd = zeros(2*NModes);

JdV = zeros(NModes,NModes*NDof);
for i = 1:NModes
    JdV(i,(i-1)*NDof+(1:NDof)) = W(:,i).'*A;
end
JdV = complex2comp(JdV,[NModes NDof]);

JdW = zeros(NModes,NModes*NDof);
for i = 1:NModes
    JdW(i,(i-1)*NDof+(1:NDof)) = (A*V(:,i)).';
end
JdW = complex2comp(JdW,[NModes NDof]);

%jacob of W.'*A*V - 1 = 0 w.r.t p
Jdp = diag(W.'*dAdp*V);
Jdp = complex2comp(Jdp,NModes);

%jacob of A*V - V*D = 0
JVd = zeros(NModes*NDof,NModes);
for i = 1:NModes
    JVd((i-1)*NDof+(1:NDof),i) = -V(:,i);
end
JVd = complex2comp(JVd,[NDof NModes]);

for i = 1:NModes
    As{i} = A - d(i)*eye(NDof);
end
JVV = blkdiag(As{:});
JVV = complex2comp(JVV,[NDof NDof]);

JVW = 0*JVV;

%jacob of A*V - V*D = 0 w.r.t p
JVp = zeros(NDof*NModes,1);
for i = 1:NModes
    JVp((i-1)*NDof+(1:NDof)) = dAdp*V(:,i);
end
JVp = complex2comp(JVp,NDof);

%jacob of W'*A - D*W' = 0
JWd = zeros(NModes*NDof,NModes);
for i = 1:NModes
    JWd((i-1)*NDof+(1:NDof),i) = -W(:,i);
end
JWd = complex2comp(JWd,[NDof NModes]);

for i = 1:NModes
    As{i} = A.' - d(i)*eye(NDof);
end
JWW = blkdiag(As{:});
JWW = complex2comp(JWW,[NDof NDof]);

JWV = 0*JWW;

%jacob of W'*A - D*W = 0 w.r.t p
JWp = zeros(NDof*NModes,1);
for i = 1:NModes
    JWp((i-1)*NDof+(1:NDof)) = dAdp.'*W(:,i);
end
JWp = complex2comp(JWp,NDof);

J = [Jdd JdV JdW Jdp;
     JVd JVV JVW JVp;
     JWd JWV JWW JWp];

J = J .* repmat(Xscale(:)',size(J,1),1);  
