function P = setupmodel(P,type)
if strcmpi(type,'rigid')
    Ar = setuprigid(P.Rotor);
elseif strcmpi(type,'FE')
    Ar = setupFE(P.Rotor);
else
    error('Unknown model type')
end

%Now add on stator DOF
As = [];
for i = 1:length(P.Stator)
    As  = [As P.Stator{i}.S'];
end
A = [Ar As];

%Impose bearing constraints
AConstr = [];
for i = 1:length(P.Bearing)
    RbSb = (P.Bearing{i}.Ro*P.Bearing{i}.So - P.Bearing{i}.Ri*P.Bearing{i}.Si);
    for k = 1:2
        if isinf(P.Bearing{i}.Kxx(k,k))
            AConstr(end+1,:) = RbSb(k,:);
        end
        if isinf(P.Bearing{i}.Kyy(k,k))
            AConstr(end+1,:) = RbSb(k+2,:);
        end
    end
end

for i = 1:length(P.Stator)
    Ss = P.Stator{i}.S;
    for k = 1:4
        if isinf(P.Stator{i}.Ks(k,k))
            AConstr(end+1,:) = Ss(k,:);
        end
    end
end

if ~isempty(AConstr)
    Af = A;
    A = Af*null(AConstr*Af,'r');
end

P.Model = enforce_constraints(P.Mesh,A);