function P = setupmodel(P,type)
if strcmpi(type,'rigid')
    RRotor = setuprigid(P.Rotor);
elseif strcmpi(type,'FE')
    RRotor = setupFE(P.Rotor,P.Bearing);
else
    error('Unknown model type')
end

%Impose bearing constraints
RBearing = [];
NBearings = length(P.Bearing);
RbBearing = zeros(2*4*NBearings,size(RRotor,1));
for i = 1:NBearings
    RbSb = (P.Bearing{i}.Ri*P.Bearing{i}.Si - P.Bearing{i}.Ro*P.Bearing{i}.So);
    RbUb = (P.Bearing{i}.Ri*P.Bearing{i}.Ui - P.Bearing{i}.Ro*P.Bearing{i}.Uo);
    for k = 1:2
        if isinf(P.Bearing{i}.Kxx(k,k))
            RBearing(end+1,:) = RbSb(k,:);
            RbBearing(:,end+1) = RbUb(k,:);
        end
        if isinf(P.Bearing{i}.Kyy(k,k))
            RBearing(end+1,:) = RbSb(k+2,:);
            RbBearing(:,end+1) = RbUb(k+2,:);
        end
    end
end

RStator = [];
for i = 1:length(P.Stator)
    Ss = P.Stator{i}.S;
    for k = 1:4
        if isinf(P.Stator{i}.Ks(k,k))
            RStator(end+1,:) = Ss(k,:);
            RbBearing(:,end+1) = 0;
        end
    end
end

AConstr = [RRotor; RBearing; RStator];
if ~isempty(AConstr)
    A = null(AConstr,'r');
else
    A = eye(P.Mesh.NDof);
end

P.Model = enforce_constraints(P.Mesh,A);

P.Mesh.R = [RRotor; RBearing; RStator]';
P.Mesh.Rotor.R    = [RRotor; 0*RBearing; 0*RStator]';
P.Mesh.Bearing.R  = [0*RRotor; RBearing; 0*RStator]';
P.Mesh.Bearing.Rb = RbBearing;
P.Mesh.Stator.R   = [0*RRotor; 0*RBearing; RStator]';
