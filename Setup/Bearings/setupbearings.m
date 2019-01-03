function B = setupbearings(B,R,O,x0)
if ~iscell(B)
    B = {B};
end
if nargin > 3
    for i = 1:length(B)
        for j = 1:2
            qi{i}{j} = B{i}.Ri{j} * B{i}.Si{j} * x0;
            qo{i}{j} = B{i}.Ro{j} * B{i}.So{j} * x0;
        end
    end
else
    for i = 1:length(B)
        for j = 1:2
            qi{i}{j} = zeros(4,1);
            qo{i}{j} = zeros(4,1);
        end
    end
end
iInputCount = 0;
iInternalCount = 0;
for i = 1:length(B)
    %if connected to ground, pad the iRotor vector
    if length(B{i}.iRotor) == 1
        B{i}.iRotor = [NaN B{i}.iRotor];
        B{i}.iNode  = [NaN B{i}.iNode];
    end
    
    if ~isnan(B{i}.iRotor(1))
        Oo = R{B{i}.iRotor(1)}.Speed*O;
    else
        Oo = 0;
    end
    
    if ~isnan(B{i}.iRotor(2))
        Oi = R{B{i}.iRotor(2)}.Speed*O;
    else
        Oi = 0;
    end
        
    B{i} = setup_each_bearing(B{i},i,R,qi{i},qo{i},Oi,Oo);
    for j = 1:2
        B{i}.iInputo{j} = iInputCount + (1:4);
        B{i}.iInputi{j} = iInputCount + (5:8);
        iInputCount = iInputCount + 8;
    end
    
    for j = 1:2
        B{i}.iInternal{j} = iInternalCount + (1:B{i}.NDofInt(j));
        iInternalCount = iInternalCount + B{i}.NDofInt(j);
    end
end

function B = setup_each_bearing(B,ind,R,xi,xo,Oi,Oo)
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

%error if we don't have the connection details
params_required = {'iRotor','iNode'};
for i = 1:length(params_required)
    if ~isfield(B,params_required{i})
        error('Cannot find parameter "%s" for %s',params_required{i},B.Name);
    end
end

%convert any linear stiffness properties (kxx) into stiffness matrices (Kxx)
dof = {'xx','xy','yy'};
mat = {'K','C'};
for iMat = 1:2
    for iDof = 1:3
        if ~isfield(B,[mat{iMat} dof{iDof}]) && isfield(B,[lower(mat{iMat}) dof{iDof}])
            if ~iscell(B.([lower(mat{iMat}) dof{iDof}]))
                B.([lower(mat{iMat}) dof{iDof}]) = {B.([lower(mat{iMat}) dof{iDof}])};
            end
            for j = 1:length(B.([lower(mat{iMat}) dof{iDof}]))
                B.([mat{iMat} dof{iDof}]){j} = diag([B.([lower(mat{iMat}) dof{iDof}]){j} 0]);
            end
        end
    end
end

if ~isfield(B,'ui')
    B.ui = {zeros(4,1),zeros(4,1)};
end
if ~isfield(B,'uo')
    B.uo = {zeros(4,1),zeros(4,1)};
end

%now extend everything to 2
for iMat = 1:2
    for iDof = 1:3
        if isfield(B,[mat{iMat} dof{iDof}])
            if ~iscell(B.([mat{iMat} dof{iDof}]))
                B.([mat{iMat} dof{iDof}]) = {B.([mat{iMat} dof{iDof}])};
            end
            if length(B.([mat{iMat} dof{iDof}])) == 1
                if iMat == 1 %k
                    B.([mat{iMat} dof{iDof}]){end+1} = diag([inf inf]);
                else %c
                    B.([mat{iMat} dof{iDof}]){end+1} = zeros(2);
                end
            end
        else
            if iDof ~= 2
                B.([mat{iMat} dof{iDof}]) = {diag([inf inf]) diag([inf inf])};
            else
                B.([mat{iMat} dof{iDof}]) = {zeros(2) zeros(2)};  
            end
        end
    end
end

model_fields = {'Model','Params'};
for iField = 1:2
    if ~iscell(B.(model_fields{iField}))
        B.(model_fields{iField}) = {B.(model_fields{iField})};
    end
    if length(B.(model_fields{iField})) == 1
        B.(model_fields{iField}){end+1} = '';
    end
end

RBear = eye(4);
RBear = RBear([1 4 2 3],:);
RBear(4,:) = -RBear(4,:);

for j = 1:2
    if strncmp(B.Model{j},'REB',3)
        [B.Params{j},B.Fb{j},B.Kb{j},B.Cb{j},B.xInt{j}] = setupREB(B.Params{j},xi{j}+B.ui{j},xo{j}+B.uo{j},Oi,Oo); 
        B.NDofInt(j) = B.Params{j}.Model.NDofTot;       

        B.Kxx{j} = B.Kb{j}(1:2,1:2);
        B.Kxy{j} = B.Kb{j}(1:2,3:4);
        B.Kyy{j} = B.Kb{j}(3:4,3:4);
        
        B.Cxx{j} = B.Cb{j}(1:2,1:2);
        B.Cxy{j} = B.Cb{j}(1:2,3:4);
        B.Cyy{j} = B.Cb{j}(3:4,3:4);
    elseif strncmp(B.Model{j},'SFD',3)
        [B.Params{j},B.Fb{j},B.Kb{j},B.Cb{j},B.xInt{j}] = setupSFD(B.Params{j},xi{j}+B.ui{j},xo{j}+B.uo{j},Oi,Oo);
        B.NDofInt(j) = 0;

        B.Kxx{j} = B.Kb{j}(1:2,1:2);
        B.Kxy{j} = B.Kb{j}(1:2,3:4);
        B.Kyy{j} = B.Kb{j}(3:4,3:4);
        
        B.Cxx{j} = B.Cb{j}(1:2,1:2);
        B.Cxy{j} = B.Cb{j}(1:2,3:4);
        B.Cyy{j} = B.Cb{j}(3:4,3:4);
    else
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
                [B.(params2default{i})] = deal(repmat({zeros(2)},1,2));
            end
        end
          
        B.NDofInt(j) = 0;
               
        B.Kb{j} = [B.Kxx{j} B.Kxy{j};
                   B.Kxy{j} B.Kyy{j}];
        
        B.Kb{j} = kron([1 -1; -1 1],B.Kb{j});
               
        B.Cb{j} = [B.Cxx{j} B.Cxy{j};
                   B.Cxy{j} B.Cyy{j}];
       
        B.Cb{j} = kron([1 -1; -1 1],B.Cb{j});
        
        Kb = max(-1E20,min(B.Kb{j},1E20));
        B.Fb{j} = Kb*[xi{j}+B.ui{j};xo{j}+B.uo{j}];
        
        B.xInt{j} = [];
    end
    B.Mb{j} = zeros(8);
    B.R{j} = blkdiag(RBear,RBear);
    B.Ri{j} = RBear;
    B.Ro{j} = RBear;
end

if isfield(B,'m') && ~isfield(B,'mx')
    B.mx = B.m;
    B.my = B.m;
end

if isfield(B,'I') && ~isfield(B,'Ixx')
    B.Ixx = B.I;
    B.Iyy = B.I;
end

params2default = {'mx','my','Ixx','Iyy'}; %mass terms
for i = 1:length(params2default)
    if ~isfield(B,params2default{i})
        [B.(params2default{i})] = deal(0);
    end
end

%setup the transformation and stiffness/damping matrices for each
%bearing
B.M = diag([B.mx B.my B.Ixx B.Iyy]);

%store the z position
B.z = R{B.iRotor(end)}.Nodes(B.iNode(end));

%throw an error if the z coords of the nodes on different rotors
%don't match
for j = 1:2
    if ~isnan(B.iRotor)
        if abs(B.z - R{B.iRotor(j)}.Nodes(B.iNode(j))) > 1E-6
            error('The z position of Rotor %d, Node %d does not match that of Rotor %d, Node %d, leading to a misalignment for Bearing %d',B.iRotor(j),B.iNode(j),B.iRotor(1),B.iNode(1),i);
        end
    end

end

%we can't have infs in the final bearing stiffness matrix as this breaks
%the FE setup code, so we set to zero now. don't worry, we'll check Kxx etc
%later to find the infs, so these values won't contribute to the final
%stiffness matrix later anyway
for j = 1:2
    ii = isinf(B.Kb{j}); B.Kb{j}(ii) = 0;
    ii = isinf(B.Cb{j}); B.Cb{j}(ii) = 0;
    %         B.Kb{j}  = max(min(B.Kb{j},1E20),-1E20);
    %         B.Cb{j}  = max(min(B.Cb{j},1E20),-1E20);
end