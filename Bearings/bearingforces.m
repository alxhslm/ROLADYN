function [Forces,Stiffness] = bearingforces(P,States)
NPts = size(States.x,2);

%default any missing fields
States = default_inputs(States,P,NPts);
States = default_speeds(States,NPts);
States = default_int(States,P,NPts);

%initialise outputs
Forces = forces_init(P,NPts);
if nargout > 1
    Stiffness = stiffness_init(P,NPts);
end

for i = 1:length(P.Bearing)
    %loop over each bearing
    
    Oshaft = repmat(0*States.O,3,1);
    Ashaft = repmat(0*States.A,3,1);
    for j = 1:2
        Oshaft(j,:) = P.Bearing{i}.Node{j}.Speed * States.O;
        Ashaft(j,:) = P.Bearing{i}.Node{j}.Speed * States.A;
    end
    
    %assemble inputs for current bearing
    StatesB = states_init_j(P.Bearing{i},Oshaft,Ashaft,States);
    
    switch P.Bearing{i}.Model
        case 'REB'
            if nargout > 1
                [ForcesB,~,StiffnessB] = REB_model(P.Bearing{i}.Params,StatesB);
                StiffnessB.K   = StiffnessB.K   + P.Bearing{i}.Params.KPar;
                StiffnessB.Kqq = StiffnessB.Kqq + P.Bearing{i}.Params.KPar;
                StiffnessB.C   = StiffnessB.C   + P.Bearing{i}.Params.CPar;
                StiffnessB.Cqq = StiffnessB.Cqq + P.Bearing{i}.Params.CPar;
            else
                ForcesB = REB_model(P.Bearing{i}.Params,StatesB);
            end
            ForcesB.F = ForcesB.F + P.Bearing{i}.Params.KPar*[StatesB.qi; StatesB.qo] + P.Bearing{i}.Params.CPar*[StatesB.qidot;StatesB.qodot];
        case 'SFD'
            if nargout > 1
                [ForcesB,~,StiffnessB] = SFD_model(P.Bearing{i}.Params,StatesB);
                StiffnessB.K = StiffnessB.K + P.Bearing{i}.Params.KSq;
            else
                ForcesB = SFD_model(P.Bearing{i}.Params,StatesB);
            end
            ForcesB.F = ForcesB.F + P.Bearing{i}.Params.KSq*[StatesB.qi; StatesB.qo];
        case 'radial'
            if nargout > 1
                [ForcesB,~,StiffnessB] = radial_model(P.Bearing{i}.Params,StatesB);
            else
                ForcesB = radial_model(P.Bearing{i}.Params,StatesB);
            end
            ForcesB.F = ForcesB.F + P.Bearing{i}.Params.KPar*[StatesB.qi; StatesB.qo];
        case 'linear'
            if nargout > 1
                [ForcesB,~,StiffnessB] = linear_model(P.Bearing{i},StatesB);
            else
                ForcesB = linear_model(P.Bearing{i},StatesB);
            end
        case 'piezo'
            if nargout > 1
                [ForcesB,~,StiffnessB] = piezo_model(P.Bearing{i}.Params,StatesB);
            else
                ForcesB = piezo_model(P.Bearing{i}.Params,StatesB);
            end
        otherwise
            if nargout > 1
                [ForcesB,~,StiffnessB] = empty_model(StatesB);
            else
                ForcesB = empty_model(StatesB);
            end
    end
    
    if ~isempty(ForcesB.FInt)
        Forces.FInt = Forces.FInt + P.Bearing{i}.V'*ForcesB.FInt;
        Forces.xInt = Forces.xInt + P.Bearing{i}.V'*ForcesB.xInt;
        Forces.xdotInt = Forces.xdotInt + P.Bearing{i}.V'*ForcesB.xdotInt;
        Forces.xddotInt = Forces.xddotInt + P.Bearing{i}.V'*ForcesB.xddotInt;
    end

    R = P.Bearing{i}.R;
    U = P.Bearing{i}.U;
    
    Forces.F  = Forces.F  + U'*R'*ForcesB.F;
    
    if nargout>1
        %assemble stiffness structure
        Stiffness.K   = Stiffness.K   + mtimesx(U'*R',mtimesx(StiffnessB.K,R*U));
        
        Stiffness.Kqq = Stiffness.Kqq + mtimesx(U'*R',mtimesx(StiffnessB.Kqq,R*U));
        Stiffness.Kqx = Stiffness.Kqx + mtimesx(U'*R',mtimesx(StiffnessB.Kqx,P.Bearing{i}.V));
        Stiffness.Kxq = Stiffness.Kxq + mtimesx(P.Bearing{i}.V',mtimesx(StiffnessB.Kxq,R*U));
        Stiffness.Kxx = Stiffness.Kxx + mtimesx(P.Bearing{i}.V',mtimesx(StiffnessB.Kxx,P.Bearing{i}.V));

        Stiffness.C   = Stiffness.C   + mtimesx(U'*R',mtimesx(StiffnessB.C,R*U));
        
        Stiffness.Cqq = Stiffness.Cqq + mtimesx(U'*R',mtimesx(StiffnessB.Cqq,R*U));
        Stiffness.Cqx = Stiffness.Cqx + mtimesx(U'*R',mtimesx(StiffnessB.Cqx,P.Bearing{i}.V));
        Stiffness.Cxq = Stiffness.Cxq + mtimesx(P.Bearing{i}.V',mtimesx(StiffnessB.Cxq,R*U));
        Stiffness.Cxx = Stiffness.Cxx + mtimesx(P.Bearing{i}.V',mtimesx(StiffnessB.Cxx,P.Bearing{i}.V));
        
        Stiffness.Mqq = Stiffness.Mqq + mtimesx(U'*R',mtimesx(StiffnessB.Mqq,R*U));
        Stiffness.Mqx = Stiffness.Mqx + mtimesx(U'*R',mtimesx(StiffnessB.Mqx,P.Bearing{i}.V));
        Stiffness.Mxq = Stiffness.Mxq + mtimesx(P.Bearing{i}.V',mtimesx(StiffnessB.Mxq,R*U));
        Stiffness.Mxx = Stiffness.Mxx + mtimesx(P.Bearing{i}.V',mtimesx(StiffnessB.Mxx,P.Bearing{i}.V));
    end
end

function StatesB = states_init_j(B,Oshaft,Ashaft,States)
StatesB.qo     = B.Ro * (B.Uo * States.x);
StatesB.qodot  = B.Ro * (B.Uo * States.xdot);
StatesB.qoddot = B.Ro * (B.Uo * States.xddot);

StatesB.qi     = B.Ri * (B.Ui * States.x);
StatesB.qidot  = B.Ri * (B.Ui * States.xdot);
StatesB.qiddot = B.Ri * (B.Ui * States.xddot);

StatesB.Oo = Oshaft(1,:); StatesB.Oi = Oshaft(2,:); 
StatesB.Ao = Ashaft(1,:); StatesB.Ai = Ashaft(2,:); 

StatesB.xInt     = B.V*States.xInt;
StatesB.xdotInt  = B.V*States.xdotInt;
StatesB.xddotInt = B.V*States.xddotInt;

StatesB.bSolve = States.bSolve;

function States = default_inputs(States,P,NPts)
if ~isfield(States,'x')
    error('Need at least the state vector x!')
end
if ~isfield(States,'xdot')
    States.xdot = 0*States.x;
end
if ~isfield(States,'xddot')
    States.xddot = 0*States.x;
end

function Forces = forces_init(P,NPts)
Forces.F = zeros(P.Mesh.Bearing.NInput,NPts);
Forces.FInt = zeros(P.Mesh.NDofInt,NPts);
Forces.xInt = zeros(P.Mesh.NDofInt,NPts);
Forces.xdotInt = zeros(P.Mesh.NDofInt,NPts);
Forces.xddotInt = zeros(P.Mesh.NDofInt,NPts);

function Stiffness = stiffness_init(P,NPts)
NBearingDof = P.Mesh.Bearing.NInput;
Stiffness.K = zeros(NBearingDof,NBearingDof,NPts);

Stiffness.Kqq = zeros(NBearingDof,NBearingDof,NPts);
Stiffness.Kqx = zeros(NBearingDof,P.Mesh.NDofInt,NPts);
Stiffness.Kxq = zeros(P.Mesh.NDofInt,NBearingDof,NPts);
Stiffness.Kxx = zeros(P.Mesh.NDofInt,P.Mesh.NDofInt,NPts);

Stiffness.C = zeros(NBearingDof,NBearingDof,NPts);

Stiffness.Cqq = zeros(NBearingDof,NBearingDof,NPts);
Stiffness.Cqx = zeros(NBearingDof,P.Mesh.NDofInt,NPts);
Stiffness.Cxq = zeros(P.Mesh.NDofInt,NBearingDof,NPts);
Stiffness.Cxx = zeros(P.Mesh.NDofInt,P.Mesh.NDofInt,NPts);

Stiffness.M = zeros(NBearingDof,NBearingDof,NPts);

Stiffness.Mqq = zeros(NBearingDof,NBearingDof,NPts);
Stiffness.Mqx = zeros(NBearingDof,P.Mesh.NDofInt,NPts);
Stiffness.Mxq = zeros(P.Mesh.NDofInt,NBearingDof,NPts);
Stiffness.Mxx = zeros(P.Mesh.NDofInt,P.Mesh.NDofInt,NPts);

function States = default_speeds(States,NPts)
fields = {'O','A'};
for i = 1:length(fields)
    if ~isfield(States,fields{i})
        States.(fields{i}) = 0;
    end
end

if length(States.O) == 1
    States.O = States.O + zeros(1,NPts);
    States.A = States.A + zeros(1,NPts);
end

function States = default_int(States,P,NPts)
if ~isfield(States,'bSolve')
    States.bSolve = 1;
end
if ~isfield(States,'xInt')
    States.xInt = zeros(P.Mesh.NDofInt,NPts);
end
if ~isfield(States,'xdotInt')
    States.xdotInt = 0*States.xInt;
end
if ~isfield(States,'xddotInt')
    States.xddotInt = 0*States.xInt;
end