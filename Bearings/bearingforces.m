function [Forces,Stiffness] = bearingforces(P,States)
NPts = size(States.x,2);

%default any missing fields
States = default_inputs(States,P,NPts);
States = default_speeds(States,NPts);
States = default_options(States);

%initialise outputs
Forces = forces_init(P,NPts);
if nargout > 1
    Stiffness = stiffness_init(P,NPts);
end

for i = 1:length(P.Bearing)
    bLin = strcmp(P.Bearing{i}.Model,'linear');
    if (States.bNL && ~bLin) ||  (States.bLin && bLin)
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
                    StiffnessB.C   = StiffnessB.C   + P.Bearing{i}.Params.CPar;
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
        
        R  = P.Bearing{i}.R;
        U  = P.Bearing{i}.U;
        Ue = P.Bearing{i}.Ue;
        
        Forces.F  = Forces.F  + U'*R'*ForcesB.F;
        
        if nargout>1
            %assemble stiffness structure
            StiffnessB = default_stiffness_j(P.Bearing{i},StiffnessB);
            
            Stiffness.K   = Stiffness.K   + mtimesx(U'*R',mtimesx(StiffnessB.K,R*U));
            Stiffness.Ku  = Stiffness.Ku  + mtimesx(U'*R',mtimesx(StiffnessB.Ku,Ue));

            Stiffness.C   = Stiffness.C   + mtimesx(U'*R',mtimesx(StiffnessB.C,R*U));   
            Stiffness.Cu  = Stiffness.Cu  + mtimesx(U'*R',mtimesx(StiffnessB.Cu,Ue));
 
            Stiffness.M   = Stiffness.M   + mtimesx(U'*R',mtimesx(StiffnessB.M,R*U));
            Stiffness.Mu  = Stiffness.Mu  + mtimesx(U'*R',mtimesx(StiffnessB.Mu,Ue));
        end
    end
end

function Stiffness = default_stiffness_j(B,Stiffness)
NPts = size(Stiffness.K,3);
NDof = 2*4;
NExc = B.NExcite;

defaultable_fields = {'C','M'};
for i = 1:2
    if ~isfield(Stiffness,defaultable_fields{i})
        Stiffness.(defaultable_fields{i}) = zeros(NDof,NDof,NPts);
    end
end

defaultable_fields = {'Ku','Cu','Mu'};
for i = 1:3
    if ~isfield(Stiffness,defaultable_fields{i})
        Stiffness.(defaultable_fields{i}) = zeros(NDof,NExc,NPts);
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

StatesB.u     = B.Ue*States.u;
StatesB.udot  = B.Ue*States.udot;
StatesB.uddot = B.Ue*States.uddot;

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
if ~isfield(States,'u')
    States.u = zeros(P.Mesh.Excite.NExcite,NPts);
end
if ~isfield(States,'udot')
    States.udot = 0*States.u;
end
if ~isfield(States,'uddot')
    States.uddot = 0*States.u;
end

function Forces = forces_init(P,NPts)
Forces.F = zeros(P.Mesh.Bearing.NInput,NPts);

function Stiffness = stiffness_init(P,NPts)
NBearingDof = P.Mesh.Bearing.NInput;
NExcite = P.Mesh.Excite.NExcite;

Stiffness.K = zeros(NBearingDof,NBearingDof,NPts);
Stiffness.Ku = zeros(NBearingDof,NExcite,NPts);

Stiffness.C = zeros(NBearingDof,NBearingDof,NPts);
Stiffness.Cu = zeros(NBearingDof,NExcite,NPts);

Stiffness.M = zeros(NBearingDof,NBearingDof,NPts);
Stiffness.Mu = zeros(NBearingDof,NExcite,NPts);

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


function States = default_options(States)
if ~isfield(States,'bNL')
    States.bNL = 1;
end
if ~isfield(States,'bLin')
    States.bLin = 1;
end