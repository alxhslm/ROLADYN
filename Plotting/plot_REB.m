function varargout = plot_REB(varargin)

nargs = nargin;
if ishandle(varargin{1})
    bPlot = 0;
    fig = varargin{1};
    varargin = varargin(2:end);
    nargs = nargs - 1;

    U = get(fig,'UserData');
    hBall = U.hBall;
    hRace = U.hRace;
else
    bPlot = 1;
end

if nargs > 5
    [REB,q,Xr,Xz,Q,Ai,Ao] = deal(varargin{:});
elseif nargs > 2
    [REB,q,Xr,Xz,Q] = deal(varargin{:});
    Ai = 0;
    Ao = 0;
else
    REB = varargin{1};
    q = zeros(4,1);
    Xz = REB.Geometry.D/2*sin(REB.Elements.alpha);
    Xr = REB.Geometry.D/2*cos(REB.Elements.alpha);
    Q = 0*Xz;
    Ai = 0;
    Ao = 0;
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
 

if bPlot
    fig = figure;
    xlabel('z (mm)');
    ylabel('x (mm)');
    zlabel('y (mm)');
    axis equal
    hold on   

    hBall = NaN(REB.Elements.N,1);
    
    if strncmpi(REB.Setup.Arrangement,'double',6)
        hRace = NaN(2,2);
    else
        hRace = NaN(2,1);
    end
end

q = R*q;
Acage = REB.Kinematics.rCagei * Ai + REB.Kinematics.rCageo * Ao;
aBall = REB.Kinematics.rPivot_lat.*(q(1).*cos(REB.Elements.psi) + q(2).*sin(REB.Elements.psi)) + REB.Kinematics.rPivot_ax*q(3) + REB.Kinematics.rPivot_rot*(q(4).*sin(REB.Elements.psi) - q(5).*cos(REB.Elements.psi)) + REB.Elements.alpha;
q(1:3,:) = q(1:3,:)*m2mm;

for i = 1:REB.Elements.N
   hBall(i) = plot_ball(hBall(i),(REB.Geometry.zo*sign(REB.Elements.z(i))+Xz(i))*m2mm,(REB.Geometry.Ro+Xr(i))*m2mm,REB.Geometry.D/2*m2mm,aBall(i),Acage+REB.Elements.psi(i),N2kN*Q(i));
end

hRace(1,1) = plot_race(hRace(1,1),  q ,REB.Geometry.zi*m2mm, REB.Geometry.ri*m2mm,REB.Geometry.Ri*m2mm,REB.Elements.z(1)*m2mm,REB.Geometry.D*m2mm,Ai,'FaceAlpha',0.5,'FaceColor','k');
hRace(2,1) = plot_race(hRace(2,1), 0*q,REB.Geometry.zo*m2mm,-REB.Geometry.ro*m2mm,REB.Geometry.Ro*m2mm,REB.Elements.z(1)*m2mm,REB.Geometry.D*m2mm,Ao,'FaceAlpha',0.5,'FaceColor','k');

if strncmpi(REB.Setup.Arrangement,'double',6)
    hRace(1,2) = plot_race(hRace(1,2),  q,-REB.Geometry.zi*m2mm, REB.Geometry.ri*m2mm,REB.Geometry.Ri*m2mm,REB.Elements.z(2)*m2mm,REB.Geometry.D*m2mm,Ai,'FaceAlpha',0.5,'FaceColor','k');
    hRace(2,2) = plot_race(hRace(2,2),0*q,-REB.Geometry.zo*m2mm,-REB.Geometry.ro*m2mm,REB.Geometry.Ro*m2mm,REB.Elements.z(2)*m2mm,REB.Geometry.D*m2mm,Ao,'FaceAlpha',0.5,'FaceColor','k');
end

U.hBall = hBall;
U.hRace = hRace;
set(fig,'UserData',U);

if bPlot
    cb = colorbar('Location','EastOutside');
    ylabel(cb,'Qi (kN)')
    caxis([0 max(Q)]*N2kN);
    view(3);
end

if nargout > 1
    varargout{1} = fig;
end

function h = plot_ball(h,z,R,r,a,psi,Q,varargin)
long = linspace(0,2*pi,10);
lat = linspace(-pi/2,pi/2,10);
[long,lat] = meshgrid(long,lat);
U = r*cos(lat).*cos(long);
V = r*sin(lat);
W = r*cos(lat).*sin(long);

%rotate by alpha
X =  U*cos(a) + V*sin(a) + z;
V = -U*sin(a) + V*cos(a) + R;

%rotate by cage position
Y = V*cos(psi) - W*sin(psi);
Z = V*sin(psi) + W*cos(psi);

C = 0*X + Q;

if isnan(h)
    h = surf(X,Y,Z,C,varargin{:});
else
    set(h,'Xdata',X,'ydata',Y,'zdata',Z,'cdata',C,varargin{:});
end

function h = plot_race(h,q,z,r,R,zB,D,A,varargin)
a1 = asin((zB + D/2 - z)/r);
a2 = asin((zB - D/2 - z)/r);
phi = linspace(a1,a2,10);
psi = linspace(0,2*pi,40) + A;
[phi,psi] = meshgrid(phi,psi);

dz = q(3)+ q(4)*R*sin(psi) - q(5)*R*cos(psi);
dr = q(1)*cos(psi) + q(2)*sin(psi) - q(4)*z*sin(psi) + q(5)*z*cos(psi);

X = z  + r*sin(phi) + dz;
Y = (R - r*cos(phi) + dr).*cos(psi);
Z = (R - r*cos(phi) + dr).*sin(psi);

C = 0*Z;

if isnan(h)
    h = surf(X,Y,Z,C,varargin{:});
else
    set(h,'Xdata',X,'ydata',Y,'zdata',Z,'cdata',C,varargin{:});
end