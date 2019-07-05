function T = enforce_constraints(S,A)
%now combine all of the rotor mass matrices into global matrices
T.M  = A'*S.M*A;
T.G  = A'*S.G*A;
T.K  = A'*S.K*A;
T.C  = A'*S.C*A;
T.Fg = A'*S.Fg;
T.F0 = A'*S.F0;
T.A = S.A*A;

%now the rotor & bearing matrices
T.Rotor.M = A'*S.Rotor.M*A;
T.Rotor.G = A'*S.Rotor.G*A;
T.Rotor.C = A'*S.Rotor.C*A;
T.Rotor.K = A'*S.Rotor.K*A;
T.Rotor.Fg = A'*S.Rotor.Fg;
T.Rotor.F0 = A'*S.Rotor.F0;

T.Bearing.K = A'*S.Bearing.K*A;
T.Bearing.C = A'*S.Bearing.C*A;
T.Bearing.M = A'*S.Bearing.M*A;
T.Bearing.Fg = A'*S.Bearing.Fg;
T.Bearing.F0 = A'*S.Bearing.F0;

%and finally the excitation matrices
T.Excite.Kgd = A'*S.Excite.Kgd;
T.Excite.Cgd = A'*S.Excite.Cgd;
T.Excite.Mgd = A'*S.Excite.Mgd;

T.Excite.Kub = A'*S.Excite.Kub;
T.Excite.Cub = A'*S.Excite.Cub;
T.Excite.Mub = A'*S.Excite.Mub;

T.Excite.Ke = A'*S.Excite.Ke;
T.Excite.Ce = A'*S.Excite.Ce;
T.Excite.Me = A'*S.Excite.Me;

%some useful numbers
T.NDof = size(T.M,1);
T.NDofInt = S.NDofInt;
T.NDofTot = T.NDof + S.NDofInt;