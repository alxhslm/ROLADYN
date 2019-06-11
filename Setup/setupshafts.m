function S = setupshafts(S,N,x0)
if ~iscell(S)
    S = {S};
end
for i = 1:length(S)
    if nargin > 2 && ~isempty(x0)
        x{i} = S{i}.S*x0;
    else
        x{i} = [];
    end
    S{i} = setup_each_shaft(S{i},N,x{i});
end

function S = setup_each_shaft(S,N,x0)
defaultable_fields = {'Element',        'Material','ri', 'ro'};
defaultable_values = {'euler_bernoulli', 'rigid',    0,    0};
for i = 1:length(defaultable_fields)
    if ~isfield(S,defaultable_fields{i})
        S.(defaultable_fields{i}) = defaultable_values{i};
    end
end

if ~isfield(S,'iNodes')
    S.iNodes = 1:length(N);
end
if ~isfield(S,'bGyro')
    S.bGyro = 1;
end

S.Material = setupmaterial(S.Material);

% params_required = {'ro'};
% for i = 1:length(params_required)
%     if ~isfield(S,params_required{i})
%         error('Cannot find parameter "%s" in the P.Rotor.Shaft structure',params_required{i});
%     end
% end

RShaft = eye(8);
RShaft = RShaft([1 4 5 8 2 3 6 7],:); %rotor to bearing coords
RShaft([6 8],:) = -1*RShaft([6 8],:); %flip signs of angles in xz plane

S.z   = N(S.iNodes);
S.Nz  = length(S.z);
S.le  = abs(diff(S.z));
if ~isempty(S.le)
    S.ze  = (S.z(1:end-1) + S.z(2:end))/2;
end

A = pi*(S.ro^2 - S.ri^2);
I = pi/4*(S.ro^4 - S.ri^4);

S.me = [];
IMap = eye(4*length(S.iNodes));
for i = 1:length(S.iNodes)
    S.SNode{i} = IMap((i-1)*4 + (1:4),:);
end
for k = 1:length(S.iNodes)-1
   [S.Ke{k},S.Me{k},S.Ge{k}] = feval(S.Element,S.Material,S.ri,S.ro,S.le(k));
   S.me(k) = S.Material.rho*A*S.le(k);
   S.R{k} = RShaft;
   if ~isempty(x0)
       S.F0{k} = S.Ke{k}*S.R{k}*S.Se{k}*x0;
   else
       S.F0{k} = zeros(8,1);
   end
   if ~S.bGyro
       S.Ge{k} = 0*S.Ge{k};
   end
   S.Se{k} = [S.SNode{k}; S.SNode{k+1}];
end

L = sum(S.le);
S.m  = sum(S.me);

S.Id = S.m/12*L^2 + S.Material.rho*I*L;
S.Ip = 2*S.Material.rho*I*L;

S.Ae = A;
S.Ie = I;

S.M = diag([S.m S.m S.Id S.Id]);
S.G = S.bGyro*blkdiag(zeros(2),antidiag(S.Ip*[1 -1]));

S.iLocal = ((S.iNodes(:)-1)*4 + (1:4))';
S.iLocal = S.iLocal(:);