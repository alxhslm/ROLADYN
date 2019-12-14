function P = setupsystem(P,type,O,A)
if ~isfield(P,'g')
    P.g = [0 0]';
end
if nargin < 2
    type = 'FE';
end
if nargin < 3
    O = 0;
end
if nargin < 4
    A = [];
end

% Rotors - first setup the mass related parameters for all of the rotors
if ~isfield(P,'Rotor')
    P.Rotor = {};
end
P.Rotor = setuprotors(P.Rotor);

% Bearings - now setup the bearing stiffness matrices etc
if ~isfield(P,'Bearing')
    P.Bearing = {};
end
P.Bearing = setupbearings(P.Bearing,P.Rotor,O);

if ~isfield(P,'Stator')
    P.Stator = {};
end
P.Stator = setupstators(P.Stator);

% Excitation - now define how the rotor will be excited
if ~isfield(P,'Excite')
    P.Excite = {};
end
P = setupexcitation(P);

%% Do the initial setup at x0 = 0

% Create mesh using specified nodes
P = setupmesh(P);

% Assemble the matrices for the rigid shaft model
P = setupmodel(P,type);

% Intial guess for x0
if isfield(P.Model,'x0') && length(P.Model.x0) == P.Model.NDof
    x0 = [P.Model.x0;
          P.Mesh.xInt];
else
    x0 = [(P.Model.K+1E5*eye(P.Model.NDof))\P.Model.Fg;
          zeros(P.Model.NDofInt,1)];
end

%% Now find equilibrium position & resetup
x0 = rotor_equib(P,x0,O,A);

P.Model.x0   = x0(1:P.Model.NDof);
P.Model.xInt = x0(P.Model.NDof+1:end);

P.Mesh.x0   = P.Model.A * P.Model.x0;
P.Mesh.xInt = P.Model.xInt;

% Bearings - now setup the bearing stiffness matrices etc
P.Rotor   = setuprotors(P.Rotor,P.Mesh.x0);
P.Bearing = setupbearings(P.Bearing,P.Rotor,O,P.Mesh.x0);

P = setuploads(P);