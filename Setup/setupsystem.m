function P = setupsystem(P,type,O)
if ~isfield(P,'g')
    P.g = [0 0]';
end
if nargin < 2
    type = 'FE';
end
if nargin < 3
    O = 0;
end

% Rotors - first setup the mass related parameters for all of the rotors
P.Rotor = setuprotors(P.Rotor);

% Bearings - now setup the bearing stiffness matrices etc
if isfield(P,'Bearing')
    P.Bearing = setupbearings(P.Bearing,P.Rotor,O);
else
    P.Bearing = {};
end

% Excitation - now define how the rotor will be excited
P = setupexcitation(P);

%% Do the initial setup at x0 = 0

% Create mesh using specified nodes
P = setupmesh(P);

% Assemble the matrices for the rigid shaft model
P = setupmodel(P,type);

% Intial guess for x0
if isfield(P.Model,'x0') && length(P.Model.x0) == P.Model.NDof
    x0 = [P.Model.x0;
        P.Mesh.Bearing.xInt];
else
    x0 = [(P.Model.K+1E5*eye(P.Model.NDof))\P.Model.Fg;
        P.Mesh.Bearing.xInt];
end

%% Now find equilibrium position & resetup
x0 = rotor_equib(P,x0,O);

P.Model.x0 = x0(1:P.Model.NDof);
P.Model.xInt = x0(P.Model.NDof+1:end);

P.Mesh.x0 = P.Model.A * x0(1:P.Model.NDof);
P.Mesh.xInt = x0(P.Model.NDof+1:end);

% Bearings - now setup the bearing stiffness matrices etc
P.Rotor = setuprotors(P.Rotor,P.Mesh.x0);
P.Bearing = setupbearings(P.Bearing,P.Rotor,O,P.Mesh.x0);

% Create mesh using specified nodes
P = setupmesh(P);

% Assemble the matrices for the rigid shaft model
P = setupmodel(P,type);