function REB = SKF1208(model,alpha)

REB.Setup.Type = 'ball';
REB.Setup.Arrangement = 'double_inline';
REB.Setup.Z = 18*2;
% REB.Setup.CbParallel = diag([0 0 0 0]);
% REB.Setup.KbParallel = diag([Inf 0 Inf 0]);

REB.Geometry.D  = 11.4E-3;
REB.Geometry.dm = 60.4E-3;
% REB.Geometry.z0 = 5.8E-3;
% REB.Geometry.alpha0 = atan(REB.Geometry.z0*2/REB.Geometry.dm);
if nargin < 2
    alpha = 10.87*pi/180;
end
    
REB.Geometry.alpha0 = alpha;
REB.Geometry.z0 = REB.Geometry.dm/2*tan(REB.Geometry.alpha0);
REB.Geometry.cr = 2.5E-6;
REB.Geometry.cz = 0E-6;

REB.Geometry.do = REB.Geometry.dm + REB.Geometry.cr + REB.Geometry.D*cos(REB.Geometry.alpha0);
REB.Geometry.di = REB.Geometry.dm - REB.Geometry.cr - REB.Geometry.D*cos(REB.Geometry.alpha0);

REB.Geometry.ro = REB.Geometry.do/cos(REB.Geometry.alpha0); %0.5*REB.dm/cos(REB.alpha0)+REB.D/2;
REB.Geometry.ri = REB.Geometry.D/2 + 1E-3;

% REB.Geometry.psi0 = (pi/2 - pi/REB.Setup.Z*floor((max(REB.Setup.Z-5,0)*0.5 + 1)));

%steel
REB.Material.v = 0.3;
REB.Material.E = 200E9;
REB.Material.rho = 8000;
 
%ISO VG32 oil
REB.Fluid.eta0 = 0.02;
REB.Fluid.alpha = 1E-8;

REB.Race.Inner.w = 18E-3;
REB.Race.Inner.t = 5E-3;

REB.Race.Outer.w = 18E-3;
REB.Race.Outer.t = 5E-3;
 
REB.Options.Control = 'outer';
REB.Options.bComplexDynLoads = 0; %doesn't work if true
REB.Options.bRaceCompliance = 0;
REB.Options.bCentrifugal = 0;
REB.Options.bGyro = 0;
REB.Options.bVC = 0; %doesn't work with c.f. loads

if nargin < 1
    model = 'harris';
end
REB.Model.Name = model;

% REB = setupREB(REB);
% plot_REB(REB)