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
        if ~isnan(P.Bearing{i}.iRotor(j))
            Oshaft(2*j-1,:) = P.Rotor{P.Bearing{i}.iRotor(j)}.Speed * States.O;
            Ashaft(2*j-1,:) = P.Rotor{P.Bearing{i}.iRotor(j)}.Speed * States.A;
        else
            Oshaft(2*j-1,:) = 0;
            Ashaft(2*j-1,:) = 0;
        end
    end
    
    for j = 1:2
       
        %assemble inputs for current bearing
        StatesB = states_init_j(P.Bearing{i},j,Oshaft,Ashaft,States);
               
        switch P.Bearing{i}.Model{j}
            case 'REB'
                if nargout > 1
                    [ForcesB,~,StiffnessB] = REB_model(P.Bearing{i}.Params{j},StatesB);
                    StiffnessB.K   = StiffnessB.K   + P.Bearing{i}.Params{j}.KPar;
                    StiffnessB.Kqq = StiffnessB.Kqq + P.Bearing{i}.Params{j}.KPar;
                    StiffnessB.C   = StiffnessB.C   + P.Bearing{i}.Params{j}.CPar;
                    StiffnessB.Cqq = StiffnessB.Cqq + P.Bearing{i}.Params{j}.CPar;
                else
                    ForcesB = REB_model(P.Bearing{i}.Params{j},StatesB);
                end
                ForcesB.F = ForcesB.F + P.Bearing{i}.Params{j}.KPar*[StatesB.qi; StatesB.qo] + P.Bearing{i}.Params{j}.CPar*[StatesB.qidot;StatesB.qodot];
            case 'SFD'
                if nargout > 1
                    [ForcesB,~,StiffnessB] = SFD_model(P.Bearing{i}.Params{j},StatesB);
                    StiffnessB.K = StiffnessB.K + P.Bearing{i}.Params{j}.KSq;
                else
                    ForcesB = SFD_model(P.Bearing{i}.Params{j},StatesB);
                end
                ForcesB.F = ForcesB.F + P.Bearing{i}.Params{j}.KSq*[StatesB.qi; StatesB.qo];
            case 'radial'
                if nargout > 1
                    [ForcesB,~,StiffnessB] = radial_model(P.Bearing{i}.Params{j},StatesB);
                else
                    ForcesB = radial_model(P.Bearing{i}.Params{j},StatesB);
                end
                ForcesB.F = ForcesB.F + P.Bearing{i}.Params{j}.KPar*[StatesB.qi; StatesB.qo];
            case 'linear'
                if nargout > 1
                    [ForcesB,~,StiffnessB] = linear_model(P.Bearing{i},j,StatesB);
                else
                    ForcesB = linear_model(P.Bearing{i},j,StatesB);
                end
            otherwise
                if nargout > 1
                    [ForcesB,~,StiffnessB] = empty_model(StatesB);
                else
                    ForcesB = empty_model(StatesB);
                end
        end
           
        if ~isempty(ForcesB.FInt)
            Forces.FInt = Forces.FInt + P.Bearing{i}.V{j}'*ForcesB.FInt;
            Forces.xInt = Forces.xInt + P.Bearing{i}.V{j}'*ForcesB.xInt;
            Forces.xdotInt = Forces.xdotInt + P.Bearing{i}.V{j}'*ForcesB.xdotInt;
        end
        
        Rq = [P.Bearing{i}.Ri{j} * P.Bearing{i}.Si{j};
              P.Bearing{i}.Ro{j} * P.Bearing{i}.So{j}];
        Ru = [P.Bearing{i}.Ri{j} * P.Bearing{i}.Ui{j};
              P.Bearing{i}.Ro{j} * P.Bearing{i}.Uo{j}];
        Forces.F  = Forces.F  + Rq'*ForcesB.F;
        Forces.Fb = Forces.Fb + Ru'*ForcesB.F;
                        
        if nargout>1
            %assemble stiffness structure      
            Stiffness.K   = Stiffness.K   + mtimesx(Rq',mtimesx(StiffnessB.K,Rq));

            Stiffness.Kqq = Stiffness.Kqq + mtimesx(Rq',mtimesx(StiffnessB.Kqq,Rq));
            Stiffness.Kqx = Stiffness.Kqx + mtimesx(Rq',mtimesx(StiffnessB.Kqx,P.Bearing{i}.V{j}));
            Stiffness.Kxq = Stiffness.Kxq + mtimesx(P.Bearing{i}.V{j}',mtimesx(StiffnessB.Kxq,Rq));
            Stiffness.Kxx = Stiffness.Kxx + mtimesx(P.Bearing{i}.V{j}',mtimesx(StiffnessB.Kxx,P.Bearing{i}.V{j}));
            
            Stiffness.Kqu = Stiffness.Kqu + mtimesx(Rq',mtimesx(StiffnessB.Kqq,Ru));
            Stiffness.Kxu = Stiffness.Kxu + mtimesx(P.Bearing{i}.V{j}',mtimesx(StiffnessB.Kxq,Ru));   
            
            Stiffness.C   = Stiffness.C   + mtimesx(Rq',mtimesx(StiffnessB.C,Rq));
            
            Stiffness.Cqq = Stiffness.Cqq + mtimesx(Rq',mtimesx(StiffnessB.Cqq,Rq));
            Stiffness.Cqx = Stiffness.Cqx + mtimesx(Rq',mtimesx(StiffnessB.Cqx,P.Bearing{i}.V{j}));
            Stiffness.Cxq = Stiffness.Cxq + mtimesx(P.Bearing{i}.V{j}',mtimesx(StiffnessB.Cxq,Rq));
            Stiffness.Cxx = Stiffness.Cxx + mtimesx(P.Bearing{i}.V{j}',mtimesx(StiffnessB.Cxx,P.Bearing{i}.V{j}));
                      
            Stiffness.Cqu = Stiffness.Cqu + mtimesx(Rq',mtimesx(StiffnessB.Cqq,Ru));
            Stiffness.Cxu = Stiffness.Cxu + mtimesx(P.Bearing{i}.V{j}',mtimesx(StiffnessB.Cxq,Ru));   
        end
    end
end

function StatesB = states_init_j(B,j,Oshaft,Ashaft,States)
StatesB.qo     = B.Ro{j} * (B.So{j} * States.x     + B.Uo{j} * States.u);
StatesB.qodot  = B.Ro{j} * (B.So{j} * States.xdot  + B.Uo{j} * States.udot);
StatesB.qoddot = B.Ro{j} * (B.So{j} * States.xddot + B.Uo{j} * States.uddot);

StatesB.qi     = B.Ri{j} * (B.Si{j} * States.x     + B.Ui{j} * States.u);
StatesB.qidot  = B.Ri{j} * (B.Si{j} * States.xdot  + B.Ui{j} * States.udot);
StatesB.qiddot = B.Ri{j} * (B.Si{j} * States.xddot + B.Ui{j} * States.uddot);

StatesB.Oo = Oshaft(j,:); StatesB.Oi = Oshaft(j+1,:); 
StatesB.Ao = Ashaft(j,:); StatesB.Ai = Ashaft(j+1,:); 

StatesB.xInt     = B.V{j}*States.xInt;
StatesB.xdotInt  = B.V{j}*States.xdotInt;
StatesB.xddotInt = B.V{j}*States.xddotInt;

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
if ~isfield(States,'u')
    States.u = zeros(size(P.Mesh.Excite.Sgd,1),NPts);
end
if ~isfield(States,'udot')
    States.udot = 0*States.u;
end
if ~isfield(States,'uddot')
    States.uddot = 0*States.u;
end

function Forces = forces_init(P,NPts)
Forces.F = zeros(P.Mesh.NDof,NPts);
Forces.Fb = zeros(2*2*4*length(P.Bearing),NPts);
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

Stiffness.Kqu = zeros(P.Mesh.NDof,2*2*4*length(P.Bearing),NPts);
Stiffness.Kxu = zeros(P.Mesh.NDofInt,2*2*4*length(P.Bearing),NPts);

Stiffness.C = zeros(P.Mesh.NDof,P.Mesh.NDof,NPts);

Stiffness.Cqq = zeros(P.Mesh.NDof,P.Mesh.NDof,NPts);
Stiffness.Cqx = zeros(P.Mesh.NDof,P.Mesh.NDofInt,NPts);
Stiffness.Cxq = zeros(P.Mesh.NDofInt,P.Mesh.NDof,NPts);
Stiffness.Cxx = zeros(P.Mesh.NDofInt,P.Mesh.NDofInt,NPts);

Stiffness.Cqu = zeros(P.Mesh.NDof,2*2*4*length(P.Bearing),NPts);
Stiffness.Cxu = zeros(P.Mesh.NDofInt,2*2*4*length(P.Bearing),NPts);

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