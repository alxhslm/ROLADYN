function [t,x,u,f,yEnd] = rotor_ode(P,ti,Oi,wi,A,y0)

if length(Oi) == 1
    Oi = 0*ti + Oi;
end

if length(wi) == 1
    wi = 0*ti + wi;
end

P.Model.bCompressREB = 0;
M = blkdiag(P.Model.M,eye(P.Model.NDof),zeros(P.Model.NDofInt),eye(2));

if nargin < 6 || isempty(y0)
    x0 = rotor_equib(P,[P.Model.x0; P.Model.xInt],Oi(1),0);
    xdot0  = 0*P.Model.x0(1:P.Model.NDof);
    y0 = [xdot0; x0;0;0];
end
ydot0 = odefun(0,y0,P,ti,Oi,wi);

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

ode_fun = @(t,y)odefun(t,y,P,ti,Oi,wi);
plot_fun = @(t,y,flag)odeplot(t,y,flag,P);
jacob_fun = @(t,y)odejacob(t,y,P,ti,Oi,wi);
options = odeset('OutputFcn',plot_fun,...
                 'Jacobian',jacob_fun,...
                 'Vectorized','on',...
                 'Mass',M,'MassSingular','maybe',...
                 'RelTol',1E-12,'AbsTol',1E-12);
             
[t,y] = ode15s(ode_fun,ti,y0,options);

% options = odeset('OutputFcn',fun,'Vectorized','on','AbsTol',1E-3);
% e0 = odefuni(0,y0,ydot0,P,O,w,U,M);
% [t,y] = ode15i(@(t,y,yp)odefuni(t,y,yp,P,O,w,U,M),time,y0,ydot0,options);

y = y'; t = t';
[ydot,Forces] = odefun(t,y,P,ti,Oi,wi);
O = ydot(end-1,:);
w = ydot(end-1,:);
A = y(end-1,:);
th = y(end-1,:);
[u,udot,uddot] = excitation_ode(P,O,w,A,th); 
f = [Forces.Fb; Forces.FInt];

[xCG,xdotCG,xInt]  = get_states(y,P);
x = [xCG; xInt]';
u = u';
f = f';

yEnd  = y(1:end-2,end);

function status = odeplot(t,y,flag,P)
persistent fig ax h peaks xdot_last count 
NDof = P.Model.NDof;
status = 0;
switch flag
    case 'init'
        fig = figure;
        [~,xdot]  = get_states(y,P);

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
        x  = get_states(y,P);
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

function e = odefuni(t,y,yp,P,O,w,U,M)
ydot = odefun(t,y,P,O,w,U);
e = M*yp - ydot;

function varargout = odefun(t,y,P,ti,Oi,wi)
[xCG,xdotCG,xInt,A,th]  = get_states(y,P);
xddotCG = 0*xCG;
x0 = P.Model.x0*(xCG(1,:)*0+1);

O = interp1(ti,Oi,t);
w = interp1(ti,wi,t);
[u,udot,uddot] = excitation_ode(P,O,w,A,th); 

%and now compute the bearing forces
States.x     = P.Model.A*xCG;
States.xdot  = P.Model.A*xdotCG;
States.xddot = P.Model.A*xddotCG;
States.xInt  = xInt;

States.A = A;
States.O = O;
States.bSolve = 0;
[Forces] = bearingforces(P,States);
Fi  = Forces.FInt;
if P.Model.bNL
    Fb  = P.Model.A'*Forces.F;
else
    Fb = P.Model.K*xCG + P.Model.C*xdotCG;
end

Fe = P.Model.Excite.M*uddot + P.Model.Excite.C*udot + P.Model.Excite.K*u;
Fr = P.Model.Rotor.K*(xCG-x0) + (P.Model.Rotor.G)*(O.*xdotCG) + P.Model.Rotor.F0;
Fg = P.Model.Fg;

%and finally compute the derivatives
Mydot = [(Fe + Fg - Fr - Fb);
            xdotCG;
            Fi;
            O;
            w];
                    
if any(isnan(Mydot) | isinf(Mydot))
    1
end
varargout{1} = Mydot;
if nargout > 1
    varargout{2} = Forces;
end

function J = odejacob(t,y,P,ti,Oi,wi)
if 1
    J = jacobian(@(x)odefun(t,x,P,ti,Oi,wi),y);
    return
end
[xCG,xdotCG,xInt,A,th] = get_states(y,P);
xddotCG = 0*xCG;
x0 = P.Model.x0*(xCG(1,:)*0+1);

O = interp1(ti,Oi,t);
w = interp1(ti,wi,t);
[u,udot,uddot] = excitation_ode(P,O,w,A,th);

%and now compute the bearing forces
States.x     = P.Model.A*xCG;
States.xdot  = P.Model.A*xdotCG;
States.xddot = P.Model.A*xddotCG;
States.xInt  = xInt;

States.A = A;
States.O = O;
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
    
function [xCG,xdotCG,xInt,A,th]  = get_states(y,P)
NDof = P.Model.NDof;
NDofInt = P.Model.NDofInt;
xdotCG = y(1:NDof,:);
xCG    = y(NDof+(1:NDof),:);
xInt   = y(2*NDof+(1:NDofInt),:);
A = y(end-1,:); %int O dt
th = y(end,:); %int w dt

function [u,udot,uddot] = get_input(U,O,A,w,th)
u     = real(U{1}*exp(1i*A)) + real(U{2}*exp(1i*th));
udot  = real(U{1}*1i*O.*exp(1i*A)) + real(U{2}*w.*exp(1i*th));
uddot = real(-U{1}*O.^2.*exp(1i*A)) + real(-U{2}*w.^2.*exp(1i*th));