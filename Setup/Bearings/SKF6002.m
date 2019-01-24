function REB = SKF6002(model,cr)

REB.Setup.Type = 'ball';
REB.Setup.Arrangement = 'single';
REB.Setup.Z = 9;

REB.Geometry.D  = 4.76E-3;
REB.Geometry.dm = 23.50E-3;
    
REB.Geometry.alpha0 = 0;
REB.Geometry.z0 = 0;

%C2 = 2.5um, C5 = 20um
if nargin < 2
    cr = 0;
end
REB.Geometry.cr = cr;
REB.Geometry.cz = 0E-6;

REB.Geometry.do = REB.Geometry.dm + REB.Geometry.cr + REB.Geometry.D;
REB.Geometry.di = REB.Geometry.dm - REB.Geometry.cr - REB.Geometry.D;

REB.Geometry.RRaceo = 1.05*REB.Geometry.D/2;
REB.Geometry.RRacei = 1.05*REB.Geometry.D/2;

% REB.Geometry.psi0 = (pi/2 - pi/REB.Setup.Z*floor((max(REB.Setup.Z-5,0)*0.5 + 1)));

%steel
REB.Material.v = 0.3;
REB.Material.E = 200E9;
REB.Material.rho = 8000;

REB.Contact.Inner.K = 1.939E10;
REB.Contact.Outer.K = 2.055E10;
 
%ISO VG32 oil
REB.Fluid.eta0 = 0.02;
REB.Fluid.alpha = 1E-8;

REB.Race.Inner.w = 9E-3;
REB.Race.Inner.t = Inf;

REB.Race.Outer.w = 9E-3;
REB.Race.Outer.t = 2E-3;
 
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