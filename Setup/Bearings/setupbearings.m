function B = setupbearings(B,R)
if ~iscell(B)
    B = {B};
end

iInputCount = 0;
iInternalCount = 0;
for i = 1:length(B)
    B{i} = setup_each_bearing(B{i},i,R);
    B{i}.iInputi = iInputCount + (1:4);
    B{i}.iInputo = iInputCount + (5:8);
    iInputCount = iInputCount + 8;
    

    B{i}.iInternal = iInternalCount + (1:B{i}.NDofInt);
    iInternalCount = iInternalCount + B{i}.NDofInt;
end

function Node = setupnodes(Node,Rotor)
%NB. order of nodes: inner, outer 
if length(Node) == 1
    %add outer connection to ground
    ground.Type = 'ground';
    Node{end+1} = ground;
end
for j = 1:length(Node)
    switch Node{j}.Type
        case 'ground'
            Node{j}.Speed = 0;
        case 'rotor'
            Node{j}.Speed = Rotor{Node{j}.iRotor}.Speed;
        case 'stator'
            Node{j}.Speed = 0;
        otherwise
            error('Unknown connection type "%s"',Node{j}.Type)
    end
end

function B = setup_each_bearing(B,ind,R)
if ~isfield(B, 'Name')
    B.Name = sprintf('Bearing %d',ind);
end

%error if no position is specified
if ~isfield(B, 'z')
    error('Cannot find parameter "z" for bearing "%s"',B.Name);
end

%error if we don't have the connection details
if ~isfield(B,'Node')
    error('Cannot find parameter "Node" for bearing "%s"',B.Name);
end

B.Node = setupnodes(B.Node,R);

%no damping as default
damping_fields = {'cxx','cxy','cyy'};
for i = 1:length(damping_fields)
    if ~isfield(B,damping_fields{i})
        B.(damping_fields{i}) = 0;
    end
end

%or inertia
inertia_fields = {'mxx','mxy','myy'};
for i = 1:length(inertia_fields)
    if ~isfield(B,inertia_fields{i})
        B.(inertia_fields{i}) = 0;
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

%convert any linear stiffness properties (kxx) into stiffness matrices (Kxx)
dof = {'xx','xy','yy'};
mat = {'K','C','M'};
for iMat = 1:3
    for iDof = 1:3
        if ~isfield(B,[mat{iMat} dof{iDof}]) && isfield(B,[lower(mat{iMat}) dof{iDof}])
            B.([mat{iMat} dof{iDof}]) = diag([B.([lower(mat{iMat}) dof{iDof}]) 0]);
        end
    end
end

switch B.Model
    case 'REB'
        [B.Params,B.Kb,B.Cb,B.Mb] = setupREB(B.Params); 
        B.F0 = zeros(8,1);
        
        B.NDofInt = B.Params.Model.NDofTot;       
        B.bActive = B.Params.bActive;  
        B.bRigid  = B.Params.bRigid;

        B.Kxx = B.Kb(1:2,1:2);
        B.Kxy = B.Kb(1:2,3:4);
        B.Kyy = B.Kb(3:4,3:4);

        B.Cxx = B.Cb(1:2,1:2);
        B.Cxy = B.Cb(1:2,3:4);
        B.Cyy = B.Cb(3:4,3:4);
        
        B.Mxx = B.Mb(1:2,1:2);
        B.Mxy = B.Mb(1:2,3:4);
        B.Myy = B.Mb(3:4,3:4);
    case 'radial'
        [B.Params,B.Kb,B.Cb,B.Mb] = setupRadial(B.Params); 
        B.F0 = zeros(8,1);
        
        B.NDofInt = 0;
        B.bActive = B.Params.bActive;   
        B.bRigid  = B.Params.bRigid;

        B.Kxx = B.Kb(1:2,1:2);
        B.Kxy = B.Kb(1:2,3:4);
        B.Kyy = B.Kb(3:4,3:4);

        B.Cxx = B.Cb(1:2,1:2);
        B.Cxy = B.Cb(1:2,3:4);
        B.Cyy = B.Cb(3:4,3:4);
        
        B.Mxx = B.Mb(1:2,1:2);
        B.Mxy = B.Mb(1:2,3:4);
        B.Myy = B.Mb(3:4,3:4);
    case 'SFD'
        [B.Params,B.Kb,B.Cb,B.Mb] = setupSFD(B.Params);
        B.F0 = zeros(8,1);
        
        B.NDofInt = 0;
        B.bActive = true(4,1);
        B.bRigid  = B.Params.bRigid;

        B.Kxx = B.Kb(1:2,1:2);
        B.Kxy = B.Kb(1:2,3:4);
        B.Kyy = B.Kb(3:4,3:4);

        B.Cxx = B.Cb(1:2,1:2);
        B.Cxy = B.Cb(1:2,3:4);
        B.Cyy = B.Cb(3:4,3:4);
        
        B.Mxx = B.Mb(1:2,1:2);
        B.Mxy = B.Mb(1:2,3:4);
        B.Myy = B.Mb(3:4,3:4);
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
        
        params2default = {'Fx','Fy'}; %off-diagonal stiffness terms
        for i = 1:length(params2default)
            if ~isfield(B,params2default{i})
                [B.(params2default{i})] = zeros(2,1);
            end
        end

        B.NDofInt = 0;

        B.Kb = [B.Kxx B.Kxy;
                B.Kxy B.Kyy];

        B.bActive = abs(diag(B.Kb)) > 0;
        B.bRigid  = isinf(diag(B.Kb));

        B.Kb = kron([1 -1; -1 1],B.Kb);

        B.Cb = [B.Cxx B.Cxy;
                B.Cxy B.Cyy];

        B.Cb = kron([1 -1; -1 1],B.Cb);
        
        B.Mb = [B.Mxx B.Mxy;
                B.Mxy B.Myy];

        B.Mb = kron([1 -1; -1 1],B.Mb);

        B.F0 = [B.Fx;B.Fy];
        B.F0 = kron([1; -1],B.F0);
end

RBear = eye(4);
RBear = RBear([1 4 2 3],:);
RBear(4,:) = -RBear(4,:);

R_fields = {'Ri','Ro'};
for i = 1:2
    if ~isfield(B,R_fields{i})
        B.(R_fields{i}) = RBear;
    end
end

%We can't have infs in the final bearing stiffness matrix as this breaks
%the FE setup code, so we set these terms to a zero. 

%Don't worry, we'll check Kxx etc later to find the Infs to determine fixed
%nodes. The values therefore won't contribute to the final stiffness matrix.

ii = isinf(B.Kb); B.Kb(ii) = 0;
ii = isinf(B.Cb); B.Cb(ii) = 0;
