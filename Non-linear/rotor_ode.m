function [t,x,u,f,yEnd] = rotor_ode(P,time,O,w,A,y0)

U = excitation_hbm(P);
if iscell(U)
    U{1} = A*U{1};
    U{2} = A*U{2};
else
    U = U * A;
end

P.Model.bCompressREB = 0;
M = blkdiag(P.Model.M,eye(P.Model.NDof),zeros(P.Model.NDofInt));

if nargin < 6 || isempty(y0)
    x0 = rotor_equib(P,[P.Model.x0; P.Model.xInt],O,0);
    xdot0  = 0*P.Model.x0(1:P.Model.NDof);
    y0 = [xdot0; x0];
end
ydot0 = odefun(0,y0,P,O,w,U);

% Jstr = [ones(P.Model.NDof) problem.sparsity(1:P.Model.NDof,:);
%         eye(P.Model.NDof) zeros(P.Model.NDof,P.Model.NDof+P.Model.NDofInt);
%         zeros(P.Model.NDofInt,P.Model.NDof)  problem.sparsity(P.Model.NDof+1:end,:)];

% ts = mean(diff(t0));
% N = floor(time(2)/t0(end));
% T = t0;
% X = x0;
% for i = 1:N
%     T = [T; T(end) + t0 + ts];
%     X = [X; x0];
% end
% t0 = T;
% x0 = X;

ode_fun = @(t,y)odefun(t,y,P,O,w,U);
plot_fun = @(t,y,flag)odeplot(t,y,flag,P);
jacob_fun = @(t,y)odejacob(t,y,P,O,w,U);
options = odeset('OutputFcn',plot_fun,...
                 'Jacobian',jacob_fun,...
                 'Vectorized','on',...
                 'Mass',M,'MassSingular','maybe',...
                 'RelTol',1E-12,'AbsTol',1E-12);
             
[t,y] = ode15s(ode_fun,time,y0,options);

% options = odeset('OutputFcn',fun,'Vectorized','on','AbsTol',1E-3);
% e0 = odefuni(0,y0,ydot0,P,O,w,U,M);
% [t,y] = ode15i(@(t,y,yp)odefuni(t,y,yp,P,O,w,U,M),time,y0,ydot0,options);

y = y'; t = t';
[ydot,Forces] = odefun(t,y,P,O,w,U);
[u,udot,uddot] = get_input(U,O,w,t);
f = [Forces.Fb; Forces.FInt];
% figure
% plot(t0,x0);
% hold on
% plot(t,x,'o');

[xCG,xdotCG,xInt]  = get_states(y,P.Model.NDof);
x = [xCG; xInt]';
u = u';
f = f';

yEnd  = y(:,end);

function status = odeplot(t,y,flag,P)
persistent fig ax h peaks xdot_last count 
NDof = P.Model.NDof;
status = 0;
switch flag
    case 'init'
        fig = figure;
        [~,xdot]  = get_states(y,NDof);

        for i = 1:NDof
            ax(i) = subplot(NDof,1,i);
            h(i) = plot(NaN,NaN);
            hold on
            xlim(t)
        end
        count = zeros(NDof,1);
        xdot_last = xdot;
        peaks = xdot_last + NaN;
    case ''           
        x  = get_states(y,NDof);
        for i = 1:NDof
            set(h(i),'xdata',[get(h(i),'xdata') t],'ydata',[get(h(i),'ydata') x(i,:)]);
            %             for j = 1:length(t)
            %                if (xdot(i,j) * xdot_last(i) < 0) && (x(i,j)>x0(i)) %detect ZC
            %                     amp = x(i,j) - x0(i);
            %                     if abs(amp - peaks(i))<1E-10
            %                         count(i) = count(i) + 1;
            %                     end
            %                     peaks(i) = amp;
            %                 end
            %                 xdot_last(i) = xdot(i,j);
            %             end
        end
        if all(count > 4)
            status = 1;
        end
        drawnow limitrate
    case 'done'
        close(fig);
end

%
% figure
% subplot(211)
% plot(t,X(:,1:2));
% subplot(212)
% semilogy(F,abs(Xfft(:,1:2)));
% xlim([0 pi/ts])
% hold on
% plot(0.33*O*[1 1],[1E-10,1],'k')
% plot(0.5*O*[1 1],[1E-10,1],'k')
% plot(O*[1 1],[1E-10,1],'k')
% plot(2*O*[1 1],[1E-10,1],'k')
% plot(3*O*[1 1],[1E-10,1],'k')
% REB = P.Bearing{1}.Params{1};
% plot(REB.rCagei*REB.Z*O*[1 1],[1E-10,1],'k')

function e = odefuni(t,y,yp,P,O,w,U,M)
ydot = odefun(t,y,P,O,w,U);
e = M*yp - ydot;

function varargout = odefun(t,y,P,O,w,U)
[xCG,xdotCG,xInt]  = get_states(y,P.Model.NDof);
xddotCG = 0*xCG;
x0 = P.Model.x0*(xCG(1,:)*0+1);

[u,udot,uddot] = get_input(U,O,w,t);

uGnd     = P.Mesh.Excite.Sgd*u;
udotGnd  = P.Mesh.Excite.Sgd*udot;
uddotGnd = P.Mesh.Excite.Sgd*uddot;

%and now compute the bearing forces
States.x     = P.Model.A*xCG;
States.xdot  = P.Model.A*xdotCG;
States.xddot = P.Model.A*xddotCG;
States.u     = uGnd;
States.udot  = udotGnd;
States.uddot = uddotGnd;
States.xInt  = xInt;

States.A = O*t;
States.O = O + 0*t;
States.bSolve = 0;
[Forces] = bearingforces(P,States);
Fi  = Forces.FInt;
if P.Model.bNL
    Fb  = P.Model.A'*Forces.F;
else
    Fb = P.Model.K*xCG + P.Model.C*xdotCG;
end

Fe = P.Model.Excite.Mub*P.Mesh.Excite.Sub*uddot + P.Model.Excite.Cub*P.Mesh.Excite.Sub*udot + P.Model.Excite.Kub*P.Mesh.Excite.Sub*u;
Fr = P.Model.Rotor.K*(xCG-x0) + (P.Model.Rotor.G*O)*xdotCG + P.Model.Rotor.F0;
Fg = P.Model.Fg;

%and finally compute the derivatives
Mydot = [(Fe + Fg - Fr - Fb);
            xdotCG;
            Fi];
                    
if any(isnan(Mydot) | isinf(Mydot))
    1
end
varargout{1} = Mydot;
if nargout > 1
    varargout{2} = Forces;
end

function J = odejacob(t,y,P,O,w,U)
if 0
    J = jacobian(@(x)odefun(t,x,P,O,w,U),y);
    return
end
[xCG,xdotCG,xInt]  = get_states(y,P.Model.NDof);
xddotCG = 0*xCG;
x0 = P.Model.x0*(xCG(1,:)*0+1);

[u,udot,uddot] = get_input(U,O,w,t);

uGnd     = P.Mesh.Excite.Sgd*u;
udotGnd  = P.Mesh.Excite.Sgd*udot;
uddotGnd = P.Mesh.Excite.Sgd*uddot;

%and now compute the bearing forces
States.x     = P.Model.A*xCG;
States.xdot  = P.Model.A*xdotCG;
States.xddot = P.Model.A*xddotCG;
States.u     = uGnd;
States.udot  = udotGnd;
States.uddot = uddotGnd;
States.xInt  = xInt;

States.A = O*t;
States.O = O + 0*t;
States.bSolve = 0;
[Forces,Stiffness] = bearingforces(P,States);

Cr = (P.Model.Rotor.C+P.Model.Rotor.G*O);
Kr = P.Model.Rotor.K;
Cb = P.Model.A'*Stiffness.Cqq*P.Model.A;
Kb = P.Model.A'*Stiffness.Kqq*P.Model.A;

%and finally compute the derivatives
J     = [  -(Cr+Cb)                     -(Kr+Kb)        -P.Model.A'*Stiffness.Kqx;
        eye(P.Model.NDof)           zeros(P.Model.NDof)    zeros(P.Model.NDof,P.Model.NDofInt);
        Stiffness.Cxq*P.Model.A    Stiffness.Kxq*P.Model.A        Stiffness.Kxx];
    
function [xCG,xdotCG,xInt]  = get_states(y,NDof)
xdotCG = y(1:NDof,:);
xCG    = y(NDof+(1:NDof),:);
xInt   = y((2*NDof+1):end,:);

function [u,udot,uddot] = get_input(U,O,w,t)
u     = real(U{1}*exp(1i*O*t)) + real(U{2}*exp(1i*w*t));
udot  = real(1i*O*U{1}*exp(1i*O*t)) + real(1i*w*U{2}*exp(1i*w*t));
uddot = real(-O^2*U{1}*exp(1i*O*t)) + real(-w^2*U{2}*exp(1i*w*t));

% u = u * (t > 0.1);
% udot = udot * (t > 0.1);
% uddot = uddot * (t > 0.1);