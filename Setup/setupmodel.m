function P = setupmodel(P,type)
if strcmpi(type,'rigid')
    Ar = setuprigid(P.Rotor);
elseif strcmpi(type,'FE')
    Ar = setupFE(P.Rotor);
else
    error('Unknown model type')
end

%Now add on bearing DOF
Ab = [];
for i = 1:length(P.Bearing)
    Ab  = [Ab P.Bearing{i}.Sb'];
end
A = [Ar Ab];

%Impose bearing constraints
AConstr = [];
for i = 1:length(P.Bearing)
    for j = 1:2
        RbSb = (P.Bearing{i}.Ro{j}*P.Bearing{i}.So{j} - P.Bearing{i}.Ri{j}*P.Bearing{i}.Si{j});
        for k = 1:2
            if isinf(P.Bearing{i}.Kxx{j}(k,k))
                AConstr(end+1,:) = RbSb(k,:);
            end
            if isinf(P.Bearing{i}.Kyy{j}(k,k))
                AConstr(end+1,:) = RbSb(k+2,:);
            end
        end
    end
end
if ~isempty(AConstr)
    Af = A;
    A = Af*null(AConstr*Af,'r');
end

%now combine all of the rotor mass matrices into global matrices
M   = A'*(P.Mesh.Rotor.M  + P.Mesh.Bearing.M)*A;
G   = A'*(P.Mesh.Rotor.G                    )*A;
K   = A'*(P.Mesh.Rotor.K  + P.Mesh.Bearing.K)*A;
C   = A'*(P.Mesh.Rotor.C  + P.Mesh.Bearing.C)*A;
Fg  = A'*(P.Mesh.Rotor.Fg + P.Mesh.Bearing.Fg);
F0  = A'*(P.Mesh.Rotor.F0 + P.Mesh.Bearing.F0);

%store the matrices
P.Model.M  = M;
P.Model.G  = G;
P.Model.K  = K;
P.Model.C  = C;
P.Model.Fg  = Fg;
P.Model.F0  = F0;
P.Model.A = A;

%now the rotor & bearing matrices
P.Model.Rotor.M = A'*P.Mesh.Rotor.M*A;
P.Model.Rotor.G = A'*P.Mesh.Rotor.G*A;
P.Model.Rotor.C = A'*P.Mesh.Rotor.C*A;
P.Model.Rotor.K = A'*P.Mesh.Rotor.K*A;
P.Model.Rotor.Fg = A'*P.Mesh.Rotor.Fg;
P.Model.Rotor.F0 = A'*P.Mesh.Rotor.F0;

P.Model.Bearing.K = A'*P.Mesh.Bearing.K*A;
P.Model.Bearing.C = A'*P.Mesh.Bearing.C*A;
P.Model.Bearing.M = A'*P.Mesh.Bearing.M*A;
P.Model.Bearing.Fg = A'*P.Mesh.Bearing.Fg;
P.Model.Bearing.F0 = A'*P.Mesh.Bearing.F0;

%and finally the excitation matrices
P.Model.Excite.Kgd = A'*P.Mesh.Excite.Kgd;
P.Model.Excite.Cgd = A'*P.Mesh.Excite.Cgd;
P.Model.Excite.Mgd = A'*P.Mesh.Excite.Mgd;

P.Model.Excite.Kub = A'*P.Mesh.Excite.Kub;
P.Model.Excite.Cub = A'*P.Mesh.Excite.Cub;
P.Model.Excite.Mub = A'*P.Mesh.Excite.Mub;

P.Model.Excite.Ke = A'*P.Mesh.Excite.Ke;
P.Model.Excite.Ce = A'*P.Mesh.Excite.Ce;
P.Model.Excite.Me = A'*P.Mesh.Excite.Me;

%some useful numbers
P.Model.NDof = size(P.Model.M,1);
P.Model.NDofInt = P.Mesh.NDofInt;
P.Model.NDofTot = P.Model.NDof + P.Mesh.NDofInt;