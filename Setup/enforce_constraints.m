function T = enforce_constraints(S,A)
%now combine all of the rotor mass matrices into global matrices
T.M  = A'*S.M*A;
T.G  = A'*S.G*A;
T.K  = A'*S.K*A;
T.C  = A'*S.C*A;
T.Fg = A'*S.Fg;
T.A = S.A*A;

%now the rotor & bearing matrices
T.Rotor.M = A'*S.Rotor.M*A;
T.Rotor.G = A'*S.Rotor.G*A;
T.Rotor.C = A'*S.Rotor.C*A;
T.Rotor.K = A'*S.Rotor.K*A;
T.Rotor.Fg = A'*S.Rotor.Fg;

T.Stator.M = A'*S.Stator.M*A;
T.Stator.C = A'*S.Stator.C*A;
T.Stator.K = A'*S.Stator.K*A;
T.Stator.Fg = A'*S.Stator.Fg;

T.Bearing.F = A'*S.Bearing.F;
T.Bearing.K = A'*S.Bearing.K*A;
T.Bearing.C = A'*S.Bearing.C*A;
T.Bearing.M = A'*S.Bearing.M*A;
T.Bearing.Lin.F = A'*S.Bearing.Lin.F;
T.Bearing.Lin.K = A'*S.Bearing.Lin.K*A;
T.Bearing.Lin.C = A'*S.Bearing.Lin.C*A;
T.Bearing.Lin.M = A'*S.Bearing.Lin.M*A;
T.Bearing.S = S.Bearing.S*A;

%and finally the excitation matrices
T.Excite.K = A'*S.Excite.K;
T.Excite.C = A'*S.Excite.C;
T.Excite.M = A'*S.Excite.M;
T.Excite.uSync = S.Excite.uSync;
T.Excite.uAsync = S.Excite.uAsync;
T.Excite.NExcite = S.Excite.NExcite;

%some useful numbers
T.NDof = size(T.M,1);