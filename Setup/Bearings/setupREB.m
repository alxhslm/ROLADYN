function [B,F,K,C,xInt] = setupREB(B,qi,qo,Oi,Oo)
if nargin < 2
    qi = zeros(4,1);
    qo = zeros(4,1);
end
if nargin < 4
    Oi = 0;
    Oo = 0;
end

if strcmpi(B.Setup.Type,'ball')
    if ~isfield(B.Geometry,'D')
        error('Need ball diameter')
    end
    B.Geometry.Dz = B.Geometry.D;
    B.Geometry.Dr = B.Geometry.D;
elseif strcmpi(B.Setup.Type,'roller_spherical')
    if ~isfield(B.Geometry,'Dz')
        error('Need roller diameter')
    end
    if ~isfield(B.Geometry,'Dr')
        error('Need roller curvature radius')
    end
elseif strcmpi(B.Setup.Type,'roller_cylindrical')
    if ~isfield(B.Geometry,'D')
        error('Need roller diameter')
    end
    if ~isfield(B.Geometry,'L')
         error('Need roller length')
    end
    if ~isfield(B.Geometry,'alphaR')
         error('Need roller-included angle')
    end
else
    error('Unknown bearing type')
end

B.Options  = setupOptions(B.Options);
B.Geometry = setupGeometry(B.Geometry);

if ~isfield(B,'Elements')
    B.Elements = struct();
end
if ~isfield(B.Elements,'N')
    if B.Options.bSjovall
        error('Need to know number of points for Sjovall integral formulation');
    else
        B.Elements.N = B.Setup.Z;
    end   
end
B.Elements.r = B.Setup.Z/B.Elements.N; %force ratio
switch B.Setup.Arrangement
    case {'single','double_alternating'}
        B.Elements.psi = B.Geometry.psi0 + (0:(B.Elements.N-1))'*(2*pi/B.Elements.N);
    case 'double_inline'
        B.Elements.psi = B.Geometry.psi0 + floor(0:0.5:(B.Elements.N/2-0.5))'*(4*pi/B.Elements.N);
end
    
B.Elements = createArrangement(B.Setup.Arrangement,B.Geometry,B.Elements);

if ~isfield(B.Setup,'KbParallel')
    B.Setup.KbParallel = zeros(4);
end

if ~isfield(B.Setup,'CbParallel')
    B.Setup.CbParallel = zeros(4);
end

B.Material = setupMaterial(B.Material);

if ~isfield(B,'Contact')
    B.Contact = struct();
end
if ~isfield(B.Contact,'Inner')
    B.Contact.Inner = struct();
end
if ~isfield(B.Contact,'Outer')
    B.Contact.Outer = struct();
end
if strcmpi(B.Setup.Type,'roller_cylindrical')
    B.Contact = setupLineContacts(B.Contact,B.Geometry,B.Material,B.Fluid);
else
    B.Contact = setupPointContacts(B.Contact,B.Geometry,B.Material,B.Fluid);
end

%Race compliance
B.Race = setupRaces(B.Race,B.Geometry,B.Setup,B.Material);
if isinf(B.Race.Inner.Kax) && B.Options.bRaceCompliancei
    warning('Disabling Race compliance for the inner race')
end
if isinf(B.Race.Outer.Kax) && B.Options.bRaceComplianceo
    warning('Disabling Race compliance for the outer race')
end
B.Race.K = 1./(B.Options.bRaceCompliancei/B.Race.Inner.K + B.Options.bRaceComplianceo/B.Race.Outer.K);

%dynamic loads
B.Dynamics.mb = 1/6*pi*B.Material.rho*B.Geometry.D^3;
B.Dynamics.Jb = 1/60*pi*B.Material.rho*B.Geometry.D^5;
B.Kinematics = setupKinematics(B.Geometry);

B.Model         = setupModel(B.Model,B.Options);
B.Model.NDofTot = B.Model.NDof * B.Elements.N;

% B = setupRaceModel(B);
% B = setupDynamicModel(B);

%find stiffness/damping in equilibrium position
psi = linspace(0,2*pi/B.Elements.N/B.Kinematics.rCagei,201);
wons = 0*psi + 1;
States.qi = qi*wons;
States.qo = qo*wons;
States.Oi = Oi*wons;
States.Oo = Oo*wons;
States.Ai = psi;
States.Ao = 0*wons;
States.bSolve = 1;
[Forces,Channels,Stiffness] = REB_model(B, States);
B.F0 = mean(Forces.F,2);
B.K0 = mean(Stiffness.K,3);
B.C0 = mean(Stiffness.C,3);
B.Channels = Channels;

B.Fi0 = mean(Forces.Fi,2);
B.Fo0 = mean(Forces.Fo,2);
B.qi0 = qi;
B.qo0 = qo;
B.x0 = mean(Forces.xInt,2);

%assemble outputs
B.KPar = kron([1 -1; -1 1], max(-1E20,min(B.Setup.KbParallel,1E20)));
K = B.K0 + kron([1 -1; -1 1], B.Setup.KbParallel);

B.CPar = kron([1 -1; -1 1], max(-1E20,min(B.Setup.CbParallel,1E20)));
C = B.C0 + kron([1 -1; -1 1], B.Setup.CbParallel);

F = B.F0 + B.KPar*[qi; qo];

xInt = B.x0;

function S = createArrangement(Arrangement,Geometry,S)
S.alpha = Geometry.alpha0 + 0*S.psi;
S.z = Geometry.z0 + 0*S.psi;

if strcmpi(Arrangement,'single')
    %done
elseif strncmpi(Arrangement,'double',6)
    if mod(length(S.psi),2)
        error('You need an even number of balls for a double arrangement')
    end
    %flip every other ball
    S.alpha(2:2:end) = -S.alpha(2:2:end);
    S.z(2:2:end) = -S.z(2:2:end);
else
    error('Unknown arrangement option')
end

function Model = setupModel(Model,Options)
Model.fun = str2func(['REB_', Model.Name]);
switch Model.Name
    case 'harris'
        Model.NDof = 0;
        if Options.bCentrifugal || Options.bGyro
            Model.NDof = Model.NDof + 2;
        end
        if Options.bRaceCompliancei
            Model.NDof = Model.NDof + 1;
        end
        if Options.bRaceComplianceo
            Model.NDof = Model.NDof + 1;
        end
    case 'ehd'
        Model.NDof = 2;
    case 'const_contact_fast'
        Model.NDof = 0;
    case 'const_contact'
        Model.NDof = 0;
        if Options.bCentrifugal || Options.bGyro
            Model.NDof = Model.NDof + 1;
        end
        if Options.bRaceCompliancei
            Model.NDof = Model.NDof + 1;
        end
        if Options.bRaceComplianceo
            Model.NDof = Model.NDof + 1;
        end
    case 'houpert'
        Model.NDof = 0;
    otherwise
        error('Unknown model option')
end

function Options = setupOptions(Options)

fields = {'bAnalyticalDeriv' , 1;
          'bSjovall'         , 0;
          'bPinned'          , 0;
          'bPivot'           , 0;
          'bGyro'            , 0;
          'bCentrifugal'     , 0;
          'bRaceCompliancei' , 0;
          'bRaceComplianceo' , 0;
          'bComplexDynLoads' , 0
          'bSaturate' , 0};

for i = 1:size(fields,1)
    if ~isfield(Options,fields{i,1})
        Options.(fields{i,1}) = fields{i,2};
    end
end

function Material = setupMaterial(Material)

if ~isfield(Material,'E')
    error('Need to know Young''s modulus')
end
if ~isfield(Material,'rho')
    error('Need to know density')
end
if ~isfield(Material,'v')
    Material.v = 0.3;
end

Material.Estar = Material.E/(1-Material.v^2)/2;

function Geometry = setupGeometry(Geometry)
default_fields = {'alpha0','psi0','z0'};
for i = 1:length(default_fields)
    if ~isfield(Geometry,default_fields{i})
        Geometry.(default_fields{i}) = 0;
    end
end

if isfield(Geometry,'cr') && isfield(Geometry,'dm')
    Geometry.do = Geometry.dm + Geometry.cr + Geometry.D*cos(Geometry.alpha0);
    Geometry.di = Geometry.dm - Geometry.cr - Geometry.D*cos(Geometry.alpha0);
elseif isfield(Geometry,'do') && isfield(Geometry,'di')
    Geometry.cr = 0.5*(Geometry.do - Geometry.di - 2*Geometry.D*cos(Geometry.alpha0));
    Geometry.dm = 0.5*(Geometry.di + Geometry.do);    
else
    error('Not eenough geometrical parameters specified')
end

%--location of race centres--
%radial location
Geometry.rRacei = Geometry.dm/2 + (Geometry.RRacei - Geometry.D/2)*cos(Geometry.alpha0);% - Geometry.cr/2;
Geometry.rRaceo = Geometry.dm/2 - (Geometry.RRaceo - Geometry.D/2)*cos(Geometry.alpha0);% + Geometry.cr/2;
%axial location
Geometry.zRacei = Geometry.z0 + (Geometry.RRacei - Geometry.D/2)*sin(Geometry.alpha0);% - Geometry.cz/2;
Geometry.zRaceo = Geometry.z0 - (Geometry.RRaceo - Geometry.D/2)*sin(Geometry.alpha0);% + Geometry.cz/2;

Geometry.A0 = Geometry.RRacei + Geometry.RRaceo - Geometry.D;
Geometry.gamma = Geometry.D/Geometry.dm;

function Kinematics = setupKinematics(Geometry)
%some kinematic properties - all derived assuming:
%   - same contact angle at inner/outer races
%   - this contact angle remains constant
Kinematics.rCagei = (1 - Geometry.gamma * cos(Geometry.alpha0))./2; 
Kinematics.rCageo = (1 + Geometry.gamma * cos(Geometry.alpha0))./2;

Kinematics.betao  = atan(sin(Geometry.alpha0)./(cos(Geometry.alpha0)+Geometry.gamma));
Kinematics.rRollo =  1./((cos(Geometry.alpha0)+tan(Kinematics.betao).*sin(Geometry.alpha0))./((1+Geometry.gamma*cos(Geometry.alpha0) + (cos(Geometry.alpha0)+tan(Kinematics.betao).*sin(Geometry.alpha0))./(1-Geometry.gamma*cos(Geometry.alpha0))).*(Geometry.gamma*cos(Kinematics.betao))));

Kinematics.betai  = atan(sin(Geometry.alpha0)./(cos(Geometry.alpha0)-Geometry.gamma));
Kinematics.rRolli = -1./((cos(Geometry.alpha0)+tan(Kinematics.betai).*sin(Geometry.alpha0))./((1+Geometry.gamma*cos(Geometry.alpha0) + (cos(Geometry.alpha0)+tan(Kinematics.betai).*sin(Geometry.alpha0))./(1-Geometry.gamma*cos(Geometry.alpha0))).*(Geometry.gamma*cos(Kinematics.betai))));

Kinematics.rPivot_rot = -((Geometry.dm-Geometry.D*cos(Geometry.alpha0))/2*cos(Geometry.alpha0) + (Geometry.z0+Geometry.D/2*sin(Geometry.alpha0))*sin(Geometry.alpha0)) / Geometry.D;
Kinematics.rPivot_ax = -cos(Geometry.alpha0)/Geometry.D;
Kinematics.rPivot_lat = sin(Geometry.alpha0)/Geometry.D;

function Race = setupRaces(Race,Geometry,Setup,Material)
theta = 2*pi/Setup.Z;
E = Material.E;
Race.Inner = setupRace(Race.Inner,theta,E,Geometry.di);
Race.Outer = setupRace(Race.Outer,theta,E,Geometry.do);

function Race = setupRace(Race,theta,E,d)
Race.I = Race.w * Race.t^3 / 12;
Race.A = Race.w * Race.t;

Race.Kax = E * Race.A * theta;
Race.Kfl = E * Race.I * theta;

if ~isfield(Race,'R')
    Race.R = d/2;
end
Race.K = Race.Kax/Race.R + Race.Kfl/Race.R^3;

function Contact = setupPointContacts(Contact,Geometry,Material,Fluid)
Rr = 1/(1/(0.5*Geometry.Dr) - 1/Geometry.RRacei);
if 0.5*Geometry.Dr > Geometry.RRacei
    error('The inner contact is non-conforming')
end
Rz = 1/(1/(0.5*Geometry.Dz) + 1/(0.5*Geometry.di/cos(Geometry.alpha0)));
if Geometry.Dz > (0.5*Geometry.di/cos(Geometry.alpha0))
    error('The inner contact is non-conforming')
end
Contact.Inner = setupHertzPointContact(Contact.Inner,Material,Rr,Rz);
% Contact.Inner.EHD = setupEHDcontact(Contact.Inner,Material,Fluid);

Rr = 1/(1/(0.5*Geometry.Dr) - 1/Geometry.RRaceo);
if 0.5*Geometry.Dr > Geometry.RRaceo
    error('The outer contact is non-conforming')
end
Rz = 1/(1/(0.5*Geometry.Dz) - 1/(0.5*Geometry.do/cos(Geometry.alpha0)));
if Geometry.Dz > (0.5*Geometry.do/cos(Geometry.alpha0))
    error('The outer contact is non-conforming')
end
Contact.Outer = setupHertzPointContact(Contact.Outer,Material,Rr,Rz);
% Contact.Outer.EHD = setupEHDcontact(Contact.Outer,Material,Fluid);

Contact.n = 1.5;
Contact.K = (Contact.Inner.K ^ (-1/Contact.n) + Contact.Outer.K ^ (-1/Contact.n)) ^ (-Contact.n);
Contact.lambda = (Contact.Outer.K / Contact.Inner.K)^(1/Contact.n);

function Contact = setupHertzPointContact(Contact,Material,Rx,Ry)
%x is radial (ie. inf for roller)
%y is axial

Contact.n = 3/2;

if ~isfield(Contact,'K')

    Sp = (1/Rx + 1/Ry);
    Fp = abs(1/Rx - 1/Ry) / Sp;
    k = solve_elliptical_contact(Fp);
    [F,E] = ellipke(1-1/k^2);

    astar = (2*k^2*E/pi)^(1/3);
    bstar = (2*E/pi/k)^(1/3);
    dstar = 2*F/pi * (pi/(2*k^2*E))^(1/3);

    Contact.Rx = Rx;
    Contact.Ry = Ry;
    Contact.R  = 1/Sp;
    
    Contact.k  = k;
    Contact.Kap = F;
    Contact.Eps = E;
    
    Contact.K  = (2*Sp*Material.Estar)/3 * (2/(Sp * dstar))^(3/2);
end

function k = solve_elliptical_contact(Fp)
if Fp == 0
	k = 1;
else
    k = nleqn(@(x)(Fp - hertz(x)),1.01);
end
    
function Fp = hertz(k)
[F,E] = ellipke(1-1/k^2);
Fp = ((k^2 + 1)*E - 2*F)/((k^2-1)*E);

function EHD = setupEHDcontact(Contact,Material,Fluid)
V = linspace(0.01,10000,30);
Q = logspace(-3,4,50);
[V,Q] = meshgrid(V,Q);

db = 0*V; C = 0*V;
for i = 1:size(V,2)
    [db(:,i), C(:,i)] = empirical_ehd_contact(Contact,Material,Fluid,V(:,i),Q(:,i));
end
%verify we can just use a modified Hertzian model at each speed
fun = @(p,x)hertz_contact(Contact.K,p(1),x+p(2)/1E6);
p0 = [1.5 0];
QiFit = 0*V;
for i = 1:size(V,2)
    [P(:,i), QiFit(:,i)] = modelfit(db(:,i),Q(:,i),fun,p0);
    p0 = P(:,i);
end

%now fit K,n,delta_0 etc as function of speed
nEHD = P(1,:)';
fun_n = @(p,x)(p(1)*(1-exp(-p(2)*x.^p(3))));
p0(1) = max(nEHD(:)) - Contact.n;
p0(2) = 1E-4;
p0(3) = 1;
p_n = modelfit(V(1,:)',nEHD,fun_n,p0);
EHD.n0  = 1.5;
EHD.Cn  = p_n(1);
EHD.kn  = p_n(2);
EHD.en  = p_n(3);

dEHD = P(2,:)'/1E6;
X = log(V(1,:)');
Y = log(dEHD*1E6);
p = polyfit(X,Y,1);
EHD.Cd = exp(p(2))/1E6;
EHD.ed = p(1); 

function [db c] = empirical_ehd_contact(Contact,Material,Fluid,vsum,Q)
% vrace = vrace*0 + 1E-10;
Estar = Material.Estar;

for j = 1:size(Q,2)
    us = 2*vsum(j) + 0.01;
    a = (3*Contact.R*Q(:,j)/(2*Estar)).^(1/3) * (2*Contact.Eps/(Contact.k*pi))^(1/3); 
    P = (2 * Estar * Contact.Rx) / (Fluid.eta0 * us);
    L = Fluid.alpha * 2 * Estar * P^(-0.25);   
    M = Q(:,j)/(2*Estar*Contact.Rx^2) * P^0.75;
    N = sqrt(Contact.Rx/Contact.Ry) * M;
   
    dbh = ((Contact.Kap/Contact.Eps)*a.^2)/(2*Contact.R);
%     dbh2 = (Q(:,j)/Contact.K).^(2/3);
%     plot(dbh,Q,dbh2,Q)
    db(:,j) = ehd_stiffness(N,L,dbh);

    Pd = (4*Contact.R*Q/us).*(Contact.Eps./Contact.Kap)./a;
    c(:,j) = ehd_damping(N,L,Pd);
end

function dbo = ehd_stiffness(M,L,dbh)
p = ((4-0.2*L)^7 + (3.5+0.1*L)^7)^(1/7);
q = -(0.6 + 0.6*(L+3)^-0.5);
D = (1 - p*M.^q);
dbo = D .* dbh;

function c = ehd_damping(M,L,P)
r =  0.98 - L/60;
s = -0.83 - L/125;
Cd = r.*M.^s;

c = P .* Cd;

function Contact = setupLineContacts(Contact,Geometry,Material,Fluid)
beta = Geometry.alpha0 - Geometry.alphaR;
Contact.Inner = setupHertzLineContact(Material,beta,Geometry.dm,Geometry.L);

beta = Geometry.alpha0 + Geometry.alphaR;
Contact.Outer = setupHertzLineContact(Material,beta,Geometry.dm,Geometry.L);

Contact.n = 1.074;
Contact.K = (Contact.Inner.K ^ (-1/Contact.n) + Contact.Outer.K ^ (-1/Contact.n)) ^ (-Contact.n);
Contact.lambda = (Contact.Outer.K / Contact.Inner.K)^(1/Contact.n);

function Contact = setupHertzLineContact(Material,beta,dm,L)
Estar = Material.Estar;
Contact.K = 0.2723 * Estar*L^0.8 * (cos(beta)/dm)^0.074;