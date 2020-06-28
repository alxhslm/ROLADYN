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
M = blkdiag(P.Model.M,eye(P.Model.NDof),zeros(P.Model.NDofInt),eye(2));

ode_fun = @(t,y)odefun(t,y,P,ti,Oi,wi);
plot_fun = @(t,y,flag)odeplot(t,y,flag,P);
jacob_fun = @(t,y)odejacob(t,y,P,ti,Oi,wi);
options = odeset('OutputFcn',plot_fun,...
                 'Jacobian',jacob_fun,...
                 'Vectorized','on',...
                 'Mass',M,'MassSingular','maybe',...
                 'RelTol',1E-12,'AbsTol',1E-12);
             
[t,y] = ode15s(ode_fun,ti,y0,options);

y = y'; t = t';
[ydot,Forces] = odefun(t,y,P,ti,Oi,wi);

%frequencies
O = ydot(end-1,:);
w = ydot(end-1,:);
A = y(end-1,:);
th = y(end-1,:);

%forces
f = [Forces.F; Forces.FInt];
[u,udot,uddot] = excitation_ode(P,O,w,A,th); 
fe = P.Model.Excite.M*uddot + P.Model.Excite.C*udot + P.Model.Excite.K*u;

%states
[xCG,xdotCG,xInt]  = get_states(y,P);
[~,xddotCG,~]  = get_states(ydot,P);

%% Create structure
ode.t = t';

ode.X     = [xCG; xInt]';
ode.Xdot  = [xdotCG; 0*xInt]';
ode.Xddot = [xddotCG; 0*xInt]';

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

function varargout = odefun(t,y,P,ti,Oi,wi)
if length(t) < size(y,2)
    t = t + 0*y(1,:);
end
[xCG,xdotCG,xInt,Theta,th]  = get_states(y,P);
xddotCG = 0*xCG;
x0 = P.Model.x0*(xCG(1,:)*0+1);

Omega = interp1(ti,Oi,t);
w = interp1(ti,wi,t);
[u,udot,uddot] = excitation_ode(P,Omega,w,Theta,th); 

%and now compute the bearing forces
States.x     = [xCG;xInt];
States.xdot  = [xdotCG;0*xInt];
States.xddot = [xddotCG;0*xInt];

States.u     = u;
States.udot  = udot;
States.uddot = uddot;

if P.Model.bNL
bearing_states = getbearingstates(States,P);
bearing_states.A = Theta;
bearing_states.O = Omega;
bearing_states.bSolve = 0;

Forces = bearingforces(P,bearing_states);
Fi  = Forces.FInt;
    Fb = P.Model.Bearing.S'*Forces.F;
else
    Fb = P.Model.Bearing.F0 + P.Model.Bearing.K*(xCG-x0) + P.Model.Bearing.C*xdotCG;
    Fi = [];
end

Fe = P.Model.Excite.M*uddot + P.Model.Excite.C*udot + P.Model.Excite.K*u;
Fr = P.Model.Rotor.K*(xCG-x0) + (P.Model.Rotor.G)*(Omega.*xdotCG) + P.Model.Rotor.C*xdotCG + P.Model.Rotor.F0;
Fs = P.Model.Stator.K*(xCG-x0) + P.Model.Stator.C*xdotCG + P.Model.Stator.F0;
Fg = P.Model.Fg;

%and finally compute the derivatives
Mydot = [(Fe + Fg - Fr - Fb - Fs);
            xdotCG;
            Fi;
            Omega;
            w];
                    
varargout{1} = Mydot;
if nargout > 1
    varargout{2} = Forces;
end

function J = odejacob(t,y,P,ti,Oi,wi)
[xCG,xdotCG,xInt,Theta,th] = get_states(y,P);
xddotCG = 0*xCG;

Omega = interp1(ti,Oi,t);
w = interp1(ti,wi,t);
[u,udot,uddot] = excitation_ode(P,Omega,w,Theta,th);

%and now compute the bearing forces
States.x     = [xCG;xInt];
States.xdot  = [xdotCG;0*xInt];
States.xddot = [xddotCG;0*xInt];

States.u     = u;
States.udot  = udot;
States.uddot = uddot;

bearing_states = getbearingstates(States,P);
bearing_states.A = Theta;
bearing_states.O = Omega;
bearing_states.bSolve = 0;

[~,Stiffness] = bearingforces(P,bearing_states);

Cr = P.Model.Rotor.C+P.Model.Rotor.G*Omega;
Kr = P.Model.Rotor.K;

Cs = P.Model.Stator.C;
Ks = P.Model.Stator.K;

Cb = P.Model.Bearing.S'*Stiffness.Cqq*P.Model.Bearing.S;
Kb = P.Model.Bearing.S'*Stiffness.Kqq*P.Model.Bearing.S;

%and finally compute the derivatives
J     = [  -(Cr+Cb+Cs)                  -(Kr+Kb+Ks)        -P.Model.Bearing.S'*Stiffness.Kqx;
        eye(P.Model.NDof)           zeros(P.Model.NDof)    zeros(P.Model.NDof,P.Model.NDofInt);
        Stiffness.Cxq*P.Model.Bearing.S    Stiffness.Kxq*P.Model.Bearing.S        Stiffness.Kxx];
    
J  = blkdiag(J,eye(2));

function [xCG,xdotCG,xInt,Theta,phPiezo]  = get_states(y,P)
NDof = P.Model.NDof;
NDofInt = P.Model.NDofInt;
xdotCG = y(1:NDof,:);
xCG    = y(NDof+(1:NDof),:);
xInt   = y(2*NDof+(1:NDofInt),:);
Theta = y(end-1,:); %int O dt
phPiezo = y(end,:); %int w dt