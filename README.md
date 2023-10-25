# ROLADYN
MATLAB suite for simulating the lateral dynamics of rotating machines.

## Example rotor-bearing system

Define a rotating machine with a flexible rotor as follows:
``` MATLAB
P.Rotor = {};
P.Bearing = {};
P.Stator = {};

% Rotor
P.Rotor{1}.Name = 'Main';
P.Rotor{1}.Speed = 1;
P.Rotor{1}.Nodes = [0,0.5,1.0]; % Nodes of the FEM

%% Flexible shaft
P.Rotor{1}.Shaft = {}
P.Rotor{1}.Shaft{1}.Name = 'Shaft';
P.Rotor{1}.Shaft{1}.iNodes = 1:3;
P.Rotor{1}.Shaft{1}.Section.ro = 0.0;
P.Rotor{1}.Shaft{1}.Section.ri = 0.1;
P.Rotor{1}.Shaft{1}.Material.Name = 'steel'; 
P.Rotor{1}.Shaft{1}.Material.E   = 210E9;
P.Rotor{1}.Shaft{1}.Material.rho = 7800;
P.Rotor{1}.Shaft{1}.Material.alpha = 0.0;
P.Rotor{1}.Shaft{1}.Material.beta  = 1e-6;
P.Rotor{1}.Shaft{1}.Options.Element = 'timoshenko';

%% Disc
P.Rotor{1}.Disc = {}
P.Rotor{1}.Disc{1}.Name = 'EndCap';
P.Rotor{1}.Disc{1}.iNode = 2;
P.Rotor{1}.Disc{1}.Type = 'Rigid';
P.Rotor{1}.Disc{1}.Options.bGyro = True;
P.Rotor{1}.Disc{1}.Material.Name = 'steel';
P.Rotor{1}.Disc{1}.Material.E = Inf; %Rigid disc
P.Rotor{1}.Disc{1}.Material.rho = 7800;
P.Rotor{1}.Disc{1}.Ring.Geometry.R = [0 0.3];
P.Rotor{1}.Disc{1}.Ring.Geometry.t = 15E-3;

% Bearings
P.Bearing{1}.Name = 'LeftBearing';
P.Bearing{1}.kyy = 1E3;
P.Bearing{1}.kxx = 1E3;
P.Bearing{1}.Node{1}.Type = 'rotor';
P.Bearing{1}.Node{1}.iNode = 1;
P.Bearing{1}.Node{1}.iRotor = 1;

P.Bearing{2}.Name = 'RightBearing';
P.Bearing{2}.kyy = 1E3;
P.Bearing{2}.kxx = 1E3;
P.Bearing{2}.Node{1}.Type = 'rotor';
P.Bearing{2}.Node{1}.iNode = 3;
P.Bearing{2}.Node{1}.iRotor = 1;

% Excitation 
P.Excite{1}.Name = 'Unbalance';
P.Excite{1}.Type = 'unbalance';
P.Excite{1}.Mode = 'sync';
P.Excite{1}.iRotor = 1;
P.Excite{1}.iDisc =  1;
P.Excite{1}.m     =  0.01;
P.Excite{1}.r     =  0.1;
P.Excite{1}.Angle =  0.0;

P = setupsystem(P)
```

You can visualise the system with:

``` MATLAB
plot_system(P)
```

## Modal analysis
To compute the natural frequencies at a range of speeds:

``` MATLAB
O = linspace(0,1E3,10) % Speed in rad/s
[V,d,W] = rotor_eig_cont(FE,O) % Using continuation to handle crossing modes
```

where $V$ and $W$ contain the left and right eigevectors, and $d$ contains the eigenvalues at each speed. 

To convert these eigenvalues/vectors into natural frequencies and mode shapes, run the following:

``` MATLAB
iNode = 2
[omega,zeta, modes] = eig2modes(V,d,W)
kappa = compute_whirl(P,modes,iNode)
```

where:
- `omega` contains the natural frequencies in rad/s
- `zeta` contains the non-dimensionalised damping ratios
- `kappa` contains the whirl direction, defined in terms of the motion at the central disc (node 2).

## Campbell diagram
To plot the Campbell diagram, run the following
``` MATLAB

plot_campbell(P,O,w,d,kappa,iNode) 
```

## Whirl
To plot the shape of the response at each resonance (ie each mode shape), run the following:
``` MATLAB
plot_whirl(P,O,modes,omega,kappa)
```

## Synchronous response
Compute the response at all frequencies with:

``` MATLAB
q = rotor_frf_modal(V,d,W,F,O)
```
And then plot the FRF using:

``` MATLAB
plot_frf(P,q,O)
```