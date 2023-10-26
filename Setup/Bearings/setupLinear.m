function [B,F0,K,C,M] = setupLinear(B)
B2default = {'Fx','Fy'}; %force at 0 position
for i = 1:length(B2default)
    if ~isfield(B,B2default{i})
        [B.(B2default{i})] = zeros(2,1);
    end
end
%no damping by default
damping_fields = {'cxx','cxy','cyy'};
for i = 1:length(damping_fields)
    if ~isfield(B,damping_fields{i})
        B.(damping_fields{i}) = 0;
    end
end

%or inertia by default
inertia_fields = {'mxx','mxy','myy'};
for i = 1:length(inertia_fields)
    if ~isfield(B,inertia_fields{i})
        B.(inertia_fields{i}) = 0;
    end
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

%throw error if we don't have stiffess
B_required = {'Kxx','Kyy','Cxx','Cyy'};
for i = 1:length(B_required)
    if ~isfield(B,B_required{i})
        error('Cannot find parameter "%s" in the B structure',B_required{i});
    end
end

B2default = {'Kxy','Cxy'}; %off-diagonal stiffness terms
for i = 1:length(B2default)
    if ~isfield(B,B2default{i})
        [B.(B2default{i})] = zeros(2);
    end
end

B2default = {'Fx','Fy'}; %off-diagonal stiffness terms
for i = 1:length(B2default)
    if ~isfield(B,B2default{i})
        [B.(B2default{i})] = zeros(2,1);
    end
end

B.F0 = [B.Fx;B.Fy];

B.K = [B.Kxx B.Kxy;
       B.Kxy B.Kyy];

B.C = [B.Cxx B.Cxy;
       B.Cxy B.Cyy];

B.M = [B.Mxx B.Mxy;
       B.Mxy B.Myy];

F0 = [B.F0;-B.F0];
K = kron([1 -1; -1 1],B.K);
C = kron([1 -1; -1 1],B.C);
M = kron([1 -1; -1 1],B.M);

B.bActive = abs(diag(B.K)) > 0;
B.bRigid  = isinf(diag(B.K));

% Clip Infs
B.K = max(-1E20,min(B.K,1E20));