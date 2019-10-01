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
    end
    
    Rq = [P.Bearing{i}.Ri * P.Bearing{i}.Si;
          P.Bearing{i}.Ro * P.Bearing{i}.So];
    Ru = [P.Bearing{i}.Ri * P.Bearing{i}.Ui;
          P.Bearing{i}.Ro * P.Bearing{i}.Uo];
    Forces.F  = Forces.F  + Rq'*ForcesB.F;
    Forces.Fb = Forces.Fb + Ru'*ForcesB.F;
    
    if nargout>1
        %assemble stiffness structure
        Stiffness.K   = Stiffness.K   + mtimesx(Rq',mtimesx(StiffnessB.K,Rq));
        
        Stiffness.Kqq = Stiffness.Kqq + mtimesx(Rq',mtimesx(StiffnessB.Kqq,Rq));
        Stiffness.Kqx = Stiffness.Kqx + mtimesx(Rq',mtimesx(StiffnessB.Kqx,P.Bearing{i}.V));
        Stiffness.Kxq = Stiffness.Kxq + mtimesx(P.Bearing{i}.V',mtimesx(StiffnessB.Kxq,Rq));
        Stiffness.Kxx = Stiffness.Kxx + mtimesx(P.Bearing{i}.V',mtimesx(StiffnessB.Kxx,P.Bearing{i}.V));

        Stiffness.C   = Stiffness.C   + mtimesx(Rq',mtimesx(StiffnessB.C,Rq));
        
        Stiffness.Cqq = Stiffness.Cqq + mtimesx(Rq',mtimesx(StiffnessB.Cqq,Rq));
        Stiffness.Cqx = Stiffness.Cqx + mtimesx(Rq',mtimesx(StiffnessB.Cqx,P.Bearing{i}.V));
        Stiffness.Cxq = Stiffness.Cxq + mtimesx(P.Bearing{i}.V',mtimesx(StiffnessB.Cxq,Rq));
        Stiffness.Cxx = Stiffness.Cxx + mtimesx(P.Bearing{i}.V',mtimesx(StiffnessB.Cxx,P.Bearing{i}.V));
    end
end

function StatesB = states_init_j(B,Oshaft,Ashaft,States)
StatesB.qo     = B.Ro * (B.So * States.x);
StatesB.qodot  = B.Ro * (B.So * States.xdot);
StatesB.qoddot = B.Ro * (B.So * States.xddot);

StatesB.qi     = B.Ri * (B.Si * States.x);
StatesB.qidot  = B.Ri * (B.Si * States.xdot);
StatesB.qiddot = B.Ri * (B.Si * States.xddot);

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
Forces.F = zeros(P.Mesh.NDof,NPts);
Forces.Fb = zeros(2*4*length(P.Bearing),NPts);
Forces.FInt = zeros(P.Mesh.NDofInt,NPts);
Forces.xInt = zeros(P.Mesh.NDofInt,NPts);
Forces.xdotInt = zeros(P.Mesh.NDofInt,NPts);
Forces.xddotInt = zeros(P.Mesh.NDofInt,NPts);

function Stiffness = stiffness_init(P,NPts)
Stiffness.K = zeros(P.Mesh.NDof,P.Mesh.NDof,NPts);

Stiffness.Kqq = zeros(P.Mesh.NDof,P.Mesh.NDof,NPts);
Stiffness.Kqx = zeros(P.Mesh.NDof,P.Mesh.NDofInt,NPts);
Stiffness.Kxq = zeros(P.Mesh.NDofInt,P.Mesh.NDof,NPts);
Stiffness.Kxx = zeros(P.Mesh.NDofInt,P.Mesh.NDofInt,NPts);

Stiffness.C = zeros(P.Mesh.NDof,P.Mesh.NDof,NPts);

Stiffness.Cqq = zeros(P.Mesh.NDof,P.Mesh.NDof,NPts);
Stiffness.Cqx = zeros(P.Mesh.NDof,P.Mesh.NDofInt,NPts);
Stiffness.Cxq = zeros(P.Mesh.NDofInt,P.Mesh.NDof,NPts);
Stiffness.Cxx = zeros(P.Mesh.NDofInt,P.Mesh.NDofInt,NPts);

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