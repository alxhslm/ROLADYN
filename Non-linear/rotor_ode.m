function ode = rotor_ode(P,ti,Oi,wi,y0)
P.Model.bNL = 1;

if length(Oi) == 1
    Oi = 0*ti + Oi;
end

if nargin < 4 || isempty(wi)
    wi = 0;
end

if length(wi) == 1
    wi = 0*ti + wi;
end

if nargin < 5 || isempty(y0)
    y0 = [0*P.Model.x0; P.Model.x0];
end

y0 = [y0;0;0];

P.Model.bCompressREB = 0;
M = blkdiag(P.Model.M,eye(2));

ode_fun = @(t,y)odefun(t,y,P,ti,Oi,wi);
plot_fun = @(t,y,flag)odeplot(t,y,flag,P);
jacob_fun = @(t,y)odejacob(t,y,P,ti,Oi,wi);
options = odeset('OutputFcn',plot_fun,...
                 'Jacobian',jacob_fun,...
                 'Vectorized','on',...
                 'Mass',M,'MassSingular','maybe',...
                 'RelTol',1E-4,'AbsTol',1E-10);
             
[t,y] = ode15s(ode_fun,ti,y0,options);

y = y'; t = t';
[ydot,Forces] = odefun(t,y,P,ti,Oi,wi);

%frequencies
O = ydot(end-1,:);
w = ydot(end-1,:);
A = y(end-1,:);
th = y(end-1,:);

%forces
f = Forces.F;
[u,udot,uddot] = excitation_ode(P,O,w,A,th); 
fe = P.Model.Excite.M*uddot + P.Model.Excite.C*udot + P.Model.Excite.K*u;

%states
[xCG,xdotCG]  = get_states(y,P);
[~,xddotCG,~]  = get_states(ydot,P);

%% Create structure
ode.t = t';

ode.X     = [xCG]';
ode.Xdot  = [xdotCG]';
ode.Xddot = [xddotCG]';

ode.U     = u';
ode.Udot  = udot';
ode.Uddot = uddot';

ode.F     = f';
ode.Fe    = fe';

ode.O = O';
ode.A = A';

ode.w = w';
ode.th = th';

function status = odeplot(t,y,flag,P)
persistent fig h count tp xp
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
        count = 0;
        xp = [];
        tp = [];
    case ''      
        count = count + 1;
        x  = get_states(y,P);
        tp = [tp t];
        xp = [xp x];
        if count > 100
            for i = 1:NDof
                set(h(i),'xdata',[get(h(i),'xdata') tp],'ydata',[get(h(i),'ydata') xp(i,:)]);
            end
            drawnow
            tp = [];
            xp = [];
            count = 0;
        end
    case 'done'
        close(fig);
end

function varargout = odefun(t,y,P,ti,Oi,wi)
if length(t) < size(y,2)
    t = t + 0*y(1,:);
end
[xCG,xdotCG,Theta,th]  = get_states(y,P);
xddotCG = 0*xCG;
x0 = P.Model.x0*(xCG(1,:)*0+1);

Omega = interp1q(ti,Oi,t')';
w = interp1q(ti,wi,t')';
[u,udot,uddot] = excitation_ode(P,Omega,w,Theta,th); 

%and now compute the bearing forces
States.x     = [xCG];
States.xdot  = [xdotCG];
States.xddot = [xddotCG];

States.u     = u;
States.udot  = udot;
States.uddot = uddot;

States.A = Theta;
States.O = Omega;

if P.Model.bNL
    Forces = bearingforces(P,States);
    Fb = P.Model.Bearing.S'*Forces.F;
else
    Fb = P.Model.Bearing.F0 + P.Model.Bearing.K*(xCG-x0) + P.Model.Bearing.C*xdotCG;
end

Fe = P.Model.Excite.M*uddot + P.Model.Excite.C*udot + P.Model.Excite.K*u;
Fr = P.Model.Rotor.K*(xCG-x0) + (P.Model.Rotor.G)*(Omega.*xdotCG) + P.Model.Rotor.C*xdotCG + P.Model.Rotor.F0;
Fs = P.Model.Stator.K*(xCG-x0) + P.Model.Stator.C*xdotCG + P.Model.Stator.F0;
Fg = P.Model.Fg;

%and finally compute the derivatives
Mydot = [(Fe + Fg - Fr - Fb - Fs);
            xdotCG;
            Omega;
            w];
                    
varargout{1} = Mydot;
if nargout > 1
    varargout{2} = Forces;
end

function J = odejacob(t,y,P,ti,Oi,wi)
[xCG,xdotCG,Theta,th] = get_states(y,P);
xddotCG = 0*xCG;

Omega = interp1(ti,Oi,t);
w = interp1(ti,wi,t);
[u,udot,uddot] = excitation_ode(P,Omega,w,Theta,th);

%and now compute the bearing forces
States.x     = [xCG];
States.xdot  = [xdotCG];
States.xddot = [xddotCG];

States.u     = u;
States.udot  = udot;
States.uddot = uddot;

Cr = P.Model.Rotor.C+P.Model.Rotor.G*Omega;
Kr = P.Model.Rotor.K;

Cs = P.Model.Stator.C;
Ks = P.Model.Stator.K;

States.A = Theta;
States.O = Omega;

if P.Model.bNL
    [~,Stiffness] = bearingforces(P,States);

    Cb = P.Model.Bearing.S'*Stiffness.C*P.Model.Bearing.S;
    Kb = P.Model.Bearing.S'*Stiffness.K*P.Model.Bearing.S;

else
    Cb = P.Model.Bearing.C;
    Kb = P.Model.Bearing.K;
end

J     = [  -(Cr+Cb+Cs)          -(Kr+Kb+Ks);
        eye(P.Model.NDof)   zeros(P.Model.NDof)];

    
J  = blkdiag(J,eye(2));

function [xCG,xdotCG,Theta,phPiezo]  = get_states(y,P)
NDof = P.Model.NDof;
xdotCG = y(1:NDof,:);
xCG    = y(NDof+(1:NDof),:);
Theta = y(end-1,:); %int O dt
phPiezo = y(end,:); %int w dt