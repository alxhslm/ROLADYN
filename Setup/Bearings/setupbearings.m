function B = setupbearings(B,R,O,x0)
if ~iscell(B)
    B = {B};
end
if nargin > 3
    for i = 1:length(B)
        qi{i} = B{i}.Ri * B{i}.Si * x0;
        qo{i} = B{i}.Ro * B{i}.So * x0;
    end
else
    for i = 1:length(B)
        qi{i} = zeros(4,1);
        qo{i} = zeros(4,1);
    end
end

iForceCount = 0;
iInternalCount = 0;
for i = 1:length(B)

    B{i} = setup_each_bearing(B{i},i,qi{i},qo{i},R,O);
    B{i}.iForceo = iForceCount + (1:4);
    B{i}.iForcei = iForceCount + (5:8);
    iForceCount = iForceCount + 8;
    

    B{i}.iInternal = iInternalCount + (1:B{i}.NDofInt);
    iInternalCount = iInternalCount + B{i}.NDofInt;
end

function Node = setupnodes(Node,Rotor)
if length(Node) == 1
    ground.Type = 'ground';
    Node = [{ground}; Node];
end
for j = 1:length(Node)
    switch Node{j}.Type
        case 'ground'
            Node{j}.Speed = 0;
        case 'rotor'
            Node{j}.Speed = Rotor{Node{j}.iRotor}.Speed;
        case 'stator'
            Node{j}.Speed = 0;
    end
end

function B = setup_each_bearing(B,ind,xi,xo,R,O)
%error if we don't have the connection details
if ~isfield(B,'Node')
    error('Cannot find parameter "%s" for %s','Node',B.Name);
end

B.Node = setupnodes(B.Node,R);
Oo = B.Node{1}.Speed*O;
Oi = B.Node{2}.Speed*O;

%no damping as default
damping_fields = {'cxx','cxy','cyy'};
for i = 1:length(damping_fields)
    if ~isfield(B,damping_fields{i})
        B.(damping_fields{i}) = 0;
    end
end

%default to linear bearing model
if ~isfield(B,'Model')
    B.Model = 'linear';
end

if ~isfield(B,'Params')
    B.Params = struct();
end

%thickness for plotting
if ~isfield(B,'t')
    B.t = 0.01;
end

if ~isfield(B, 'Name')
    B.Name = sprintf('Bearing %d',ind);
end

%convert any linear stiffness properties (kxx) into stiffness matrices (Kxx)
dof = {'xx','xy','yy'};
mat = {'K','C'};
for iMat = 1:2
    for iDof = 1:3
        if ~isfield(B,[mat{iMat} dof{iDof}]) && isfield(B,[lower(mat{iMat}) dof{iDof}])
            B.([mat{iMat} dof{iDof}]) = diag([B.([lower(mat{iMat}) dof{iDof}]) 0]);
        end
    end
end

if ~isfield(B,'ui')
    B.ui = zeros(4,1);
end
if ~isfield(B,'uo')
    B.uo = zeros(4,1);
end

RBear = eye(4);
RBear = RBear([1 4 2 3],:);
RBear(4,:) = -RBear(4,:);

switch B.Model
    case 'REB'
        [B.Params,B.Fb,B.Kb,B.Cb,B.xInt] = setupREB(B.Params,xi+B.ui,xo+B.uo,Oi,Oo); 
        B.NDofInt = B.Params.Model.NDofTot;       
        B.bActive = B.Params.bActive;       

        B.Kxx = B.Kb(1:2,1:2);
        B.Kxy = B.Kb(1:2,3:4);
        B.Kyy = B.Kb(3:4,3:4);

        B.Cxx = B.Cb(1:2,1:2);
        B.Cxy = B.Cb(1:2,3:4);
        B.Cyy = B.Cb(3:4,3:4);
    case 'radial'
        [B.Params,B.Fb,B.Kb,B.Cb,B.xInt] = setupRadial(B.Params,xi+B.ui,xo+B.uo,Oi,Oo); 
        B.NDofInt = 0;
        B.bActive = B.Params.bActive;   

        B.Kxx = B.Kb(1:2,1:2);
        B.Kxy = B.Kb(1:2,3:4);
        B.Kyy = B.Kb(3:4,3:4);

        B.Cxx = B.Cb(1:2,1:2);
        B.Cxy = B.Cb(1:2,3:4);
        B.Cyy = B.Cb(3:4,3:4);
    case 'SFD'
        [B.Params,B.Fb,B.Kb,B.Cb,B.xInt] = setupSFD(B.Params,xi+B.ui,xo+B.uo,Oi,Oo);
        B.NDofInt = 0;
        B.bActive = true(4,1);

        B.Kxx = B.Kb(1:2,1:2);
        B.Kxy = B.Kb(1:2,3:4);
        B.Kyy = B.Kb(3:4,3:4);

        B.Cxx = B.Cb(1:2,1:2);
        B.Cxy = B.Cb(1:2,3:4);
        B.Cyy = B.Cb(3:4,3:4);
    otherwise
        %throw error if we don't have stiffess
        params_required = {'Kxx','Kyy','Cxx','Cyy'};
        for i = 1:length(params_required)
            if ~isfield(B,params_required{i})
                error('Cannot find parameter "%s" in the B structure',params_required{i});
            end
        end

        params2default = {'Kxy','Cxy'}; %off-diagonal stiffness terms
        for i = 1:length(params2default)
            if ~isfield(B,params2default{i})
                [B.(params2default{i})] = zeros(2);
            end
        end

        B.NDofInt = 0;

        B.Kb = [B.Kxx B.Kxy;
                B.Kxy B.Kyy];

        B.bActive = abs(diag(B.Kb)) > 0;

        B.Kb = kron([1 -1; -1 1],B.Kb);

        B.Cb = [B.Cxx B.Cxy;
                   B.Cxy B.Cyy];

        B.Cb = kron([1 -1; -1 1],B.Cb);

        Kb = max(-1E20,min(B.Kb,1E20));
        B.Fb = Kb*[xi+B.ui;xo+B.uo];

        B.xInt = [];
end

B.Ri = RBear;
B.Ro = RBear;
B.R  = blkdiag(RBear,RBear);

%We can't have infs in the final bearing stiffness matrix as this breaks
%the FE setup code, so we set these terms to a small value. It must be 
%non-zero so we can detect which direction they can exert forces, for the 
%C-B model reduciton on the rotor

%Don't worry, we'll check Kxx etc later to find the Infs to determine fixed
%nodes. The values therefore won't contribute to the final stiffness matrix.

ii = isinf(B.Kb); B.Kb(ii) = 100;
ii = isinf(B.Cb); B.Cb(ii) = 100;
%         B.Kb  = max(min(B.Kb,1E20),-1E20);
%         B.Cb  = max(min(B.Cb,1E20),-1E20);
