function P = setupsystem(P,type,O,A)

%% Default inputs

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

%% Setup components

% Rotors - first setup the rotors including shaft FE models etc
if ~isfield(P,'Rotor')
    P.Rotor = {};
end
P.Rotor = setuprotors(P.Rotor);

% Bearings - now setup the bearing stiffness matrices etc
if ~isfield(P,'Bearing')
    P.Bearing = {};
end
P.Bearing = setupbearings(P.Bearing,P.Rotor);

% Stators - any non-rotating components such as bearing housings etc
if ~isfield(P,'Stator')
    P.Stator = {};
end
P.Stator = setupstators(P.Stator);

% Excitation - now define how the system will be excited
if ~isfield(P,'Excite')
    P.Excite = {};
end
P = setupexcitation(P);

%% Assembly

% Create mesh using specified nodes
P = setupmesh(P);

% Assemble the matrices and perform model reduction
P = setupmodel(P,type);

%% Find equilibrium position
P = rotor_equib(P,O,A);

%compute the loads at setup
P = setuploadandstiffness(P,O,A);