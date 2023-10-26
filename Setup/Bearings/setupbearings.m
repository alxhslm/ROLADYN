function B = setupbearings(B,R)
if ~iscell(B)
    B = {B};
end

iInputCount = 0;
for i = 1:length(B)
    B{i} = setup_each_bearing(B{i},i,R);
    B{i}.iInputi = iInputCount + (1:4);
    B{i}.iInputo = iInputCount + (5:8);
    iInputCount = iInputCount + 8;
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

if ~isfield(B,'Params')
    B.Params = struct();
end

%thickness for plotting
if ~isfield(B,'t')
    B.t = 0.01;
end

if ~isfield(B,'Model')
    error(['Bearing ' B.Name ' has no model type'])
end
switch B.Model
    case 'REB'
        B.SetupFun = str2func('setupREB');
        B.ModelFun = str2func('REB_model');
    case 'SFD'
        B.SetupFun = str2func('setupSFD');
        B.ModelFun = str2func('SFD_model');
    case 'radial'
        B.SetupFun = str2func('setupRadial');
        B.ModelFun = str2func('radial_model');
    case 'linear'
        B.SetupFun = str2func('setupLinear');
        B.ModelFun = str2func('linear_model');
    case 'empty'
        B.SetupFun = str2func('setupEmpty');
        B.ModelFun = str2func('empty_model');
    otherwise
        error(['Bearing ' B.Name ' has an invalid model type'])
end
[B.Params,B.F0,B.Kb,B.Cb,B.Mb] = B.SetupFun(B.Params); 

B.Kxx = B.Kb(1:2,1:2);
B.Kxy = B.Kb(1:2,3:4);
B.Kyy = B.Kb(3:4,3:4);

B.Cxx = B.Cb(1:2,1:2);
B.Cxy = B.Cb(1:2,3:4);
B.Cyy = B.Cb(3:4,3:4);

B.Mxx = B.Mb(1:2,1:2);
B.Mxy = B.Mb(1:2,3:4);
B.Myy = B.Mb(3:4,3:4);

B.bActive = B.Params.bActive;
B.bRigid  = B.Params.bRigid; 
        
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
    