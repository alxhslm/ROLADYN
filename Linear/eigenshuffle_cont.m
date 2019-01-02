function [V,d,W,p] = eigenshuffle_cont(fun,p,options)
p0 = p(1);
pEnd = p(end);

[A,B] = feval(fun,p0);
[V,d,W] = eig(A,B,'vector');
[V,W] = mkorthog(B,V,W);

N = size(A,1);
x0 = packeig(V,d,W);

Xscale = 0*x0 + norm(x0);
Xscale(end+1) = mean([p0,pEnd]);

pMax = max(pEnd,p0);
pMin = min(pEnd,p0);

fun_constr = @(x)eigen_cont(x,fun,Xscale);
fun_jacobian = @(x)jacobian(@eigen_cont,x,fun,Xscale);

p = p0;
x = x0;
Xprev = [x0; p0]./Xscale;
F = fun_constr(Xprev);

options = defaultmissingoptions(options);
step = options.step;
max_step = options.dpMax;
min_step = options.dpMin;
count = 0;

hWaitbar = waitbar(0, 'Frequency Range');
A = full(feval(fun_jacobian,Xprev));
tangent_prev = mynull(A);
tangent_prev = tangent_prev * sign(tangent_prev(end)) * sign(pEnd - p0);

fCamp = figure;
h = plot(p0,imag(d(1:2:end))/2/pi,'o-');
xlim([pMin wMax]);
xlabel('Speed (rad/s)')
ylabel('Freq (Hz)')
drawnow

while p(end) <= wMax && p(end) >= pMin
    X = Xprev + step*tangent_prev;
    iter = 0;
    Xlast = X + Inf;
    tangent = tangent_prev;
    while (norm(X - Xlast) > options.xtol || norm(F) > options.ftol)
        Xlast = X;
        A = full(feval(fun_jacobian,X));
        F = feval(fun_constr,X);
        X = Xlast - A\F;
        tangent = mynull(A);
 
        iter = iter + 1;
        if iter > options.maxit
             break;
        end
    end
    
    pCurr = X(end).*Xscale(end);
    
    if iter < options.maxit && (pCurr - p(end))*(pEnd - p0) > 0
        step = norm(X - Xlast);
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
            [V(:,:,end+1),d(:,end+1),W(:,:,end+1)] = unpackeig(x(:,end),N);
            
            for i = 1:N/2
                set(h(i),'xdata',p,'ydata',imag(d(2*i-1,:))/2/pi);
            end
            
        end
    else
        %failed
        step = step*options.c;
        count = count + 1;
        if count > 4
            break
        end
    end
    
    %update our progress
    waitbar((pCurr-p0)/(pEnd - p0),hWaitbar,sprintf('Step Size: %f',step));
    drawnow
end
close(hWaitbar)
close(fCamp)

function t = mynull(A)
[~,~,V] = svd(A,0);
t = V(:,end);

function x = packeig(V,d,W)
x = [V(:);d(:);W(:)];
x = [real(x); imag(x)];

function [V,d,W] = unpackeig(x,N)
x = reshape(x,[],2);
x = x(:,1) + 1i*x(:,2);
V = reshape(x(1:N*N),N,N);
d = x(N*N + (1:N));
W = reshape(x(N*N + N + (1:N*N)),N,N);

function f = eigen_cont(X,fun,Xscale)
x = X(1:end-1).*Xscale(1:end-1);
w0 = X(end).*Xscale(end);

[A,B] = feval(fun,w0);
N = size(A,1);

[V,d,W] = unpackeig(x,N);
f = eigen_gen(V,d,W,A,B);

f = [real(f);imag(f)];

function f = eigen_gen(V,d,W,A,B)
D = diag(d);
residV = (B*V*D - A*V);
residW = (D*W'*B - W'*A);
resid_norm =  diag(W'*B*V) - 1;

f = [residV(:);residW(:);resid_norm];

function [V,W] = mkorthog(B,V,W)
 %ensure we normalise the modes correctly
scale = sqrt(diag(W' * B * V)).';    
scale(scale == 0) = 1;

%now scale the eigenvectors
p = size(V,1);
V = V./ repmat(scale,p,1);
W = W./ repmat(conj(scale),p,1);

function opt = defaultmissingoptions(opt)
f = {'step','C', 'c','maxit','xtol','ftol','dwMin','dwMax'};
d = {1,     1.1, 0.5, 30    ,   1E-4,  1E-4,    0   ,  Inf};
for i = 1:length(f)
    if ~isfield(opt,f{i})
        opt.(f{i}) = d{i};
    end
end