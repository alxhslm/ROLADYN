function Model = mesh2model(Mesh,A)

%now combine all of the rotor mass matrices into global matrices
M   = A'*(Mesh.Rotor.M  + Mesh.Bearing.M)*A;
G   = A'*(Mesh.Rotor.G                    )*A;
K   = A'*(Mesh.Rotor.K  + Mesh.Bearing.K)*A;
C   = A'*(Mesh.Rotor.C  + Mesh.Bearing.C)*A;
Fg  = A'*(Mesh.Rotor.Fg + Mesh.Bearing.Fg);
F0  = A'*(Mesh.Rotor.F0 + Mesh.Bearing.F0);

%store the matrices
Model.M  = M;
Model.G  = G;
Model.K  = K;
Model.C  = C;
Model.Fg  = Fg;
Model.F0  = F0;
Model.A = A;

%now the rotor & bearing matrices
Model.Rotor.M = A'*Mesh.Rotor.M*A;
Model.Rotor.G = A'*Mesh.Rotor.G*A;
Model.Rotor.C = A'*Mesh.Rotor.C*A;
Model.Rotor.K = A'*Mesh.Rotor.K*A;
Model.Rotor.Fg = A'*Mesh.Rotor.Fg;
Model.Rotor.F0 = A'*Mesh.Rotor.F0;

Model.Bearing.K = A'*Mesh.Bearing.K*A;
Model.Bearing.C = A'*Mesh.Bearing.C*A;
Model.Bearing.M = A'*Mesh.Bearing.M*A;
Model.Bearing.Fg = A'*Mesh.Bearing.Fg;
Model.Bearing.F0 = A'*Mesh.Bearing.F0;

%and finally the excitation matrices
Model.Excite.Kgd = A'*Mesh.Excite.Kgd;
Model.Excite.Cgd = A'*Mesh.Excite.Cgd;
Model.Excite.Mgd = A'*Mesh.Excite.Mgd;

Model.Excite.Kub = A'*Mesh.Excite.Kub;
Model.Excite.Cub = A'*Mesh.Excite.Cub;
Model.Excite.Mub = A'*Mesh.Excite.Mub;

Model.Excite.Ke = A'*Mesh.Excite.Ke;
Model.Excite.Ce = A'*Mesh.Excite.Ce;
Model.Excite.Me = A'*Mesh.Excite.Me;

%some useful numbers
Model.NDof = size(Model.M,1);
Model.NDofInt = Mesh.NDofInt;
Model.NDofTot = Model.NDof + Mesh.NDofInt;