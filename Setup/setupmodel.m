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

P.Model = mesh2model(P.Mesh,A);