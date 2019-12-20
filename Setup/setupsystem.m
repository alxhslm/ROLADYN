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
P.Bearing = setupbearings(P.Bearing,P.Rotor);

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

%% Now find equilibrium position & resetup
P = rotor_equib(P,O,A);

%compute the loads at setup
P = setuploadandstiffness(P,O,A);