function fig = plot_REB_slice(varargin)

nargs = nargin;
if ishandle(varargin{1})
    bPlot = 0;
    fig = deal(varargin{1});
    varargin = varargin(2:end);
    nargs = nargs - 1;
    U = get(fig,'UserData');
else
    bPlot = 1;
    U = struct();
end

if nargs > 5
    [REB,q,U.Xr,U.Xz,U.Q,iPlot] = deal(varargin{:});
elseif nargs > 4
    [REB,q,U.Xr,U.Xz,U.Q] = deal(varargin{:});
    iPlot = 1;
else
    REB = varargin{1};
    q = zeros(4,1);
    U.Xz = REB.z;
    U.Xr = 0*REB.psi;
    U.Q = 0*U.Xz;
    iPlot = 1;
end

m2mm = 1000;
N2kN = 1E-3;

if size(q,1) == 4
    R = [1     0     0     0
         0     0     1     0
         0     0     0     0
         0     0     0    -1
         0     1     0     0
         0     0     0     0];
else
    R = eye(5);
end
 
q = R*q;
U.q = q;
U.aBall = REB.rPivot_lat.*(q(1).*cos(REB.psi) + q(2).*sin(REB.psi)) + REB.rPivot_ax*q(3) + REB.rPivot_rot*(q(4).*sin(REB.psi) - q(5).*cos(REB.psi));
U.q(1:3,:) = U.q(1:3,:)*m2mm;

if bPlot
    fig = figure;
    xlabel('z (mm)');
    ylabel('r (mm)');
    axis equal
    hold on   

    U.hBall = NaN(REB.Z,1);
    U.hRace = NaN(2,1);

    U.hControl = uicontrol('parent',fig,'Style','popupmenu','Units','Normalized','Position',[0.78  0.94  0.14  0.05],'String',num2cellstr(1:REB.Z), 'Callback', @controlcb,'Value',iPlot);       
    [U.hBall,U.hRace] = plot_ball_and_races(U.hBall,U.hRace,REB,U.q,U.Xr(iPlot),U.Xz(iPlot),U.aBall(iPlot),U.Q(iPlot),REB.z(iPlot),REB.psi(iPlot),m2mm,N2kN);
    
    U.m2mm = m2mm;
    U.N2kN = m2mm;
else
    iPlot = get(U.hControl,'Value');
    plot_ball_and_races(U.hBall,U.hRace,REB,U.q,U.Xr(iPlot),U.Xz(iPlot),U.aBall(iPlot),U.Q(iPlot),REB.z(iPlot),REB.psi(iPlot),m2mm,N2kN);
end
U.REB = REB;
set(fig,'UserData',U);

if bPlot
    cb = colorbar('Location','EastOutside');
    ylabel(cb,'Qi (kN)')
    caxis([0 max(U.Q)]*N2kN)
end

function controlcb(source,event)
i = get(source,'Value');
fig = get(source,'Parent');
U = get(fig,'UserData');
plot_ball_and_races(U.hBall,U.hRace,U.REB,U.q,U.Xr(i),U.Xz(i),U.aBall(i),U.Q(i),U.REB.z(i),U.REB.psi(i),U.m2mm,U.N2kN);

function [hBall,hRace] = plot_ball_and_races(hBall,hRace,REB,q,Xr,Xz,aBall,Q,z,psi,m2mm,N2kN)
hBall = plot_ball(hBall,(REB.zo*sign(z)+Xz)*m2mm,(REB.Ro+Xr)*m2mm,REB.D/2*m2mm,aBall,N2kN*Q);
hRace(1) = plot_race(hRace(1), q  ,sign(z)*REB.zi*m2mm, REB.ri*m2mm,REB.Ri*m2mm,z*m2mm,REB.D*m2mm,psi,'Color','k');
hRace(2) = plot_race(hRace(2), 0*q,sign(z)*REB.zo*m2mm,-REB.ro*m2mm,REB.Ro*m2mm,z*m2mm,REB.D*m2mm,psi,'Color','k');

function h = plot_ball(h,z,R,r,a,Q,varargin)
lat = linspace(0,2*pi,100); lat(end+1) = pi;
U = r*cos(lat);
V = r*sin(lat);

X =  U*cos(a) + V*sin(a) + z;
Y = -U*sin(a) + V*cos(a) + R;

C = 0*X + Q;

if isnan(h)
    h = fill(X,Y,C,varargin{:});
    h = double(h);
else
    set(h,'Xdata',X,'ydata',Y,'cdata',C,varargin{:});
end

function h = plot_race(h,q,z,r,R,zB,D,psi,varargin)
a1 = asin((zB + D/2 - z)/r);
a2 = asin((zB - D/2 - z)/r);
phi = linspace(a1,a2,100);

dz = q(3)+ q(4)*R*sin(psi) - q(5)*R*cos(psi);
dr = q(1)*cos(psi) + q(2)*sin(psi) - q(4)*z*sin(psi) + q(5)*z*cos(psi);

X = z + r*sin(phi) + dz;
Y = R - r*cos(phi) + dr;

if isnan(h)
    h = plot(X,Y,varargin{:});
    h = double(h);
else
    set(h,'Xdata',X,'ydata',Y,varargin{:});
end