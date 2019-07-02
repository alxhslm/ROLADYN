function [Forces,Channels,Stiffness] = REB_model(Params, States)
%convert bearing displacememnts into the q vector
if size(States.qi,1) == 4
    R = [1     0     0     0
         0     0     1     0
         0     0     0     0
         0     0     0    -1
         0     1     0     0
         0     0     0     0];
else
    R =  [1     0     0     0    0
          0     1     0     0    0
          0     0     1     0    0
          0     0     0     1    0
          0     0     0     0    1
          0     0     0     0    0];
end

%solve for the ball equilibrium
model = Params.Model.fun;
NPts = size(States.qi,2);

States = default_inputs(States,R);
States = default_speeds(States,NPts);
States = default_int(States,Params,NPts);

%decide if we need to show the waitbar
if ~isfield(States,'bWaitbar')
    States.bWaitbar = 0;
end

if States.bSolve
    %QS solution
    States.xInt     = zeros(Params.Model.NDofTot,NPts);
    %(Params.Geometry.ro-Params.Geometry.D/2) * repmat([sin(Params.alpha); cos(Params.alpha)],1,NPts);
    States.xdotInt  = 0*States.xInt;
    States.xddotInt = 0*States.xInt;
    States = equilibrium(model,Params,States);
end

if nargout<3
    [Forces,Channels] = feval(model,Params,States);
else
    [Forces,Channels,Stiffness] = stiffnessAndDamping(model,Params,States,R);
end
Forces.Fi = R'*Forces.Fi;
Forces.Fo = R'*Forces.Fo;
Forces.F  = [Forces.Fi;
             Forces.Fo];
             
Forces.xInt = States.xInt;
Forces.xdotInt = States.xdotInt;
Forces.xddotInt = States.xddotInt;

Forces.q = R'*(States.qi - States.qo);

function States = default_int(States,REB,NPts)
if ~isfield(States,'bSolve')
    States.bSolve = 1;
end
if ~isfield(States,'xInt')
    States.xInt = zeros(REB.Model.NDofTot,NPts);
end
if ~isfield(States,'xdotInt')
    States.xdotInt = 0*States.xInt;
end
if ~isfield(States,'xddotInt')
    States.xddotInt = 0*States.xInt;
end

function States = default_inputs(States,R)
if ~isfield(States,'qo')
    States.qo = 0*States.qi;
end
suffix = {'dot','ddot'};
IO = {'i','o'};
for i = 1:length(suffix)
    for k = 1:2
        States.(['q' IO{k} suffix{i}]) = 0*States.qi;
    end
end

prefix = {'q','q','q'};
suffix = {'' ,'dot','ddot'};
IO = {'i','o'};
for i = 1:length(prefix)
    for k = 1:2
        States.([prefix{i} IO{k} suffix{i}]) = R*States.([prefix{i} IO{k} suffix{i}]);
    end
end

function States = default_speeds(States,NPts)
fields = {'Oi','Oo','Ai','Ao'};
for i = 1:length(fields)
    if ~isfield(States,fields{i})
        States.(fields{i}) = 0;
    end
end

if length(States.Oi) == 1
    States.Oi = States.Oi + zeros(1,NPts);
    States.Oo = States.Oo + zeros(1,NPts);
    States.Ai = States.Ai + zeros(1,NPts);
    States.Ao = States.Ao + zeros(1,NPts);
end

function States = equilibrium(model,Params,States)
if Params.Model.NDof > 0
    State_j = getjthpoint(States,1);
    x0 = State_j.xInt;
        
    if States.bWaitbar
        h = waitbar(0,'Sweep');
    end
    Npts = size(States.qi,2);
    for j = 1:Npts  
        
        State_j = getjthpoint(States,j);
        
    %     Jstr = sparse(repmat(eye(B.Z) ,B.NDof,B.NDof));
        Jstr = sparse(repmat(eye(Params.Elements.N),Params.Model.NDof,Params.Model.NDof));

        funcs.objective = @(x)(0);
        funcs.gradient = @(x)(0*x0);
        funcs.constraints = @(x)equilibrium_constr(x,model,Params,State_j);
        funcs.jacobian = @(x)equilibrium_jac(x,model,Params,State_j);
        funcs.jacobianstructure = @(x)Jstr;
        funcs.iterfunc = @iteration;

        options.cl = 0*x0;
        options.cu = 0*x0;
        options.ipopt.hessian_approximation = 'limited-memory';
        options.ipopt.print_level = 0;
        options.ipopt.max_iter = 50;
        options.ipopt.constr_viol_tol = 1E-10;

        info.status = 2;
        while ~any(info.status == [0 1])
            [x, info] = ipopt(x0,funcs,options);      
            x0 = x + 1E-6*rand(size(x0)); 
        end
        States.xInt(:,j) = x;
        x0 = x;
        
        if States.bWaitbar
            waitbar(j/Npts,h);
        end
    end
    
    if States.bWaitbar
        close(h);
    end
else
    States.xInt = zeros(0,size(States.qi,2));
end
States.xdotInt = 0*States.xInt;
States.xddotInt = 0*States.xInt;

function FInt = equilibrium_constr(x,model,Params,States)
States.xInt = x;
Forces = feval(model,Params,States);
FInt = Forces.FInt;

function Kxx = equilibrium_jac(x,model,Params,States)
States.xInt = x;
[F,~,S] = feval(model,Params,States);
f0 = F.FInt;
if ~Params.Options.bAnalyticalDeriv || isempty(S)
    h = 1E-10;
    Kxx = zeros(Params.Model.NDofTot,Params.Model.NDofTot);
    x0 = x;
    for i = 1:Params.Model.NDofTot
        x = x0;
        x(i,:) = x(i,:) + h;
        f = equilibrium_constr(x,model,Params,States);
        Kxx(:,i,:) = (f-f0)./h;
    end
else
    Kxx = S.Kxx;
end
Kxx = sparse(Kxx);

function States = getjthpoint(States,j)
prefix = {'q','q','q','O','A'};
suffix = {'' ,'dot','ddot','',''};
IO = {'i','o'};
for i = 1:length(prefix)
    for k = 1:2
        States.([prefix{i} IO{k} suffix{i}]) = States.([prefix{i} IO{k} suffix{i}])(:,j);
    end
end

fields = {'xInt','xdotInt','xddotInt'};
for i = 1:length(fields)
    States.(fields{i}) = States.(fields{i})(:,j);
end

function [F,V,S] = stiffnessAndDamping(model,Params,States,R)
[F,V,S] = feval(model,Params,States);

N = size(R,2);
NPts = size(States.qi,2);
if ~Params.Options.bAnalyticalDeriv || isempty(S)
    S = struct();
    Fi0 = F.Fi;
    Fo0 = F.Fo;
    Fx = F.FInt;
    h = 1E-10;
    H = zeros(N,NPts);
    
    %%%%%%%%%%%%%%%% STIFFNESS %%%%%%%%%%%%%%%%%

    %qi derivatives
    qi = States.qi;
    S.Kqiqi = zeros(N,N,NPts);
    S.Kqoqi = zeros(N,N,NPts);
    S.Kxqi = zeros(Params.Model.NDofTot,N,NPts);  
    for i = 1:N
        H(i,:) = H(i,:) + h;
        States.qi = States.qi + R*H;
        Forces = feval(model,Params,States);
        S.Kqiqi(:,i,:) = permute(R'*(Forces.Fi-Fi0)./h,[1 3 2]);
        S.Kqoqi(:,i,:) = permute(R'*(Forces.Fo-Fo0)./h,[1 3 2]);
        S.Kxqi(:,i,:)  = permute((Forces.FInt-Fx)./h,[1 3 2]);
        States.qi = qi;
        H = 0*H;
    end
    
    %qo derivatives
    qo = States.qo;
    S.Kqiqo = zeros(N,N,NPts);
    S.Kqoqo = zeros(N,N,NPts);
    S.Kxqo = zeros(Params.Model.NDofTot,N,NPts);  
    for i = 1:N
        H(i,:) = H(i,:) + h;
        States.qo = States.qo + R*H;
        Forces = feval(model,Params,States);
        S.Kqiqo(:,i,:) = permute(R'*(Forces.Fi-Fi0)./h,[1 3 2]);
        S.Kqoqo(:,i,:) = permute(R'*(Forces.Fo-Fo0)./h,[1 3 2]);
        S.Kxqo(:,i,:) = permute((Forces.FInt-Fx)./h,[1 3 2]);
        States.qo = qo;
        H = 0*H;
    end
    
    %xInt derivatives
    x = States.xInt;
    S.Kqix = zeros(N,Params.Model.NDofTot,NPts);
    S.Kqox = zeros(N,Params.Model.NDofTot,NPts);
    S.Kxx = zeros(Params.Model.NDofTot,Params.Model.NDofTot,NPts);
    for i = 1:size(x,1)
        States.xInt(i,:) = States.xInt(i,:) + h;
        Forces = feval(model,Params,States);
        S.Kqix(:,i,:) = permute(R'*(Forces.Fi-Fi0)./h,[1 3 2]);
        S.Kqox(:,i,:) = permute(R'*(Forces.Fo-Fo0)./h,[1 3 2]);
        S.Kxx(:,i,:) = permute((Forces.FInt-Fx)./h,[1 3 2]);
        States.xInt = x;
    end
        
    %%%%%%%%%%%%%%%% DAMPING %%%%%%%%%%%%%%%%%
    
    %qi derivatives
    qidot = States.qidot;
    S.Cqiqi = zeros(N,N,NPts);
    S.Cqoqi = zeros(N,N,NPts);
    S.Cxqi = zeros(Params.Model.NDofTot,N,NPts);  
    for i = 1:N
        H(i,:) = H(i,:) + h;
        States.qidot = States.qidot + R*H;
        Forces = feval(model,Params,States);
        S.Cqiqi(:,i,:) = permute(R'*(Forces.Fi-Fi0)./h,[1 3 2]);
        S.Cqoqi(:,i,:) = permute(R'*(Forces.Fo-Fo0)./h,[1 3 2]);
        S.Cxqi(:,i,:) = permute((Forces.FInt-Fx)./h,[1 3 2]);
        States.qidot = qidot;
        H = 0*H;
    end
    
    %qo derivatives
    qodot = States.qodot;
    S.Cqiqo = zeros(N,N,NPts);
    S.Cqoqo = zeros(N,N,NPts);
    S.Cxqo = zeros(Params.Model.NDofTot,N,NPts);  
    for i = 1:N
        H(i,:) = H(i,:) + h;
        States.qodot = States.qodot + R*H;
        Forces = feval(model,Params,States);
        S.Cqiqo(:,i,:) = permute(R'*(Forces.Fi-Fi0)./h,[1 3 2]);
        S.Cqoqo(:,i,:) = permute(R'*(Forces.Fo-Fo0)./h,[1 3 2]);
        S.Cxqi(:,i,:) = permute((Forces.FInt-Fx)./h,[1 3 2]);
        States.qodot = qodot;
        H = 0*H;
    end
    
    %xInt derivatives
    xdot = States.xdotInt;
    S.Cqix = zeros(N,Params.Model.NDofTot,NPts);
    S.Cqox = zeros(N,Params.Model.NDofTot,NPts);
    S.Cxx = zeros(Params.Model.NDofTot,Params.Model.NDofTot,NPts);
    States.xdotInt = xdot;
    for i = 1:size(xdot,1)
        States.xdotInt(i,:) = States.xdotInt(i,:) + h;
        Forces = feval(model,Params,States);
        S.Cqix(:,i,:) = permute(R'*(Forces.Fi-Fi0)./h,[1 3 2]);
        S.Cqox(:,i,:) = permute(R'*(Forces.Fo-Fo0)./h,[1 3 2]);
        S.Cxx(:,i,:) = permute((Forces.FInt-Fx)./h,[1 3 2]);
        States.xdotInt = xdot;
    end
    
    
else
    %just need to multiply by R to make sure the sizes are consistent
    S.Kqiqi = mtimesx(R',mtimesx(S.Kqiqi,R));
    S.Kqoqi = mtimesx(R',mtimesx(S.Kqoqi,R));
    S.Kqiqo = mtimesx(R',mtimesx(S.Kqiqo,R));
    S.Kqoqo = mtimesx(R',mtimesx(S.Kqoqo,R));
    
    S.Kxqi = mtimesx(S.Kxqi,R);
    S.Kxqo = mtimesx(S.Kxqo,R);
    
    S.Kqix = mtimesx(R',S.Kqix);
    S.Kqox = mtimesx(R',S.Kqox);
           
    S.Cqiqi = mtimesx(R',mtimesx(S.Cqiqi,R));
    S.Cqoqi = mtimesx(R',mtimesx(S.Cqoqi,R));
    S.Cqiqo = mtimesx(R',mtimesx(S.Cqiqo,R));
    S.Cqoqo = mtimesx(R',mtimesx(S.Cqoqo,R));
    
    S.Cxqi = mtimesx(S.Cxqi,R);
    S.Cxqo = mtimesx(S.Cxqo,R);
    
    S.Cqix = mtimesx(R',S.Cqix);
    S.Cqox = mtimesx(R',S.Cqox);
end

%Combine stiffness terms
S.Kqq = [S.Kqiqi S.Kqiqo; 
         S.Kqoqi S.Kqoqo];
S.Kqx = [S.Kqix;
         S.Kqox];
S.Kxq = [S.Kxqi S.Kxqo];

%now work out effective stiffness
Kxxinv = ball_inv(S.Kxx,Params.Model.NDof);
if Params.Model.NDof > 0
    S.K = S.Kqq - mtimesx(S.Kqx,mtimesx(Kxxinv,S.Kxq));
else
    S.K = S.Kqq;
end

S.Kii = S.K(1:N,1:N,:);
S.Kio = S.K(1:N,N+(1:N),:);
S.Koi = S.K(N+(1:N),1:N,:);
S.Koo = S.K(N+(1:N),N+(1:N),:);

S.Cqq = [S.Cqiqi S.Cqiqo; 
         S.Cqoqi S.Cqoqo];
S.Cqx = [S.Cqix; 
         S.Cqox];
S.Cxq = [S.Cxqi S.Cxqo];

%now work out effective damping
Cxxinv = ball_inv(S.Cxx,Params.Model.NDof);
if Params.Model.NDof > 0
    S.C = S.Cqq - mtimesx(S.Cqx,mtimesx(Cxxinv,S.Cxq));
else
    S.C = S.Cqq;
end

S.Cii = S.C(1:N,1:N,:);
S.Cio = S.C(1:N,N+(1:N),:);
S.Coi = S.C(N+(1:N),1:N,:);
S.Coo = S.C(N+(1:N),N+(1:N),:);

function J = myinv(M)
if size(M,1) == 1
    J = 1./M;
    J(M == 0) = 0;
    return;
end
a = M(1,1,:);
b = M(1,2,:);
c = M(2,1,:);
d = M(2,2,:);

det = (a.*d - b.*c) + eps;

J = [d./det -b./det;
    -c./det  a./det];
          
function J = ball_inv(M,NDof)
if isempty(M)
    J = M;
    return;
end
Nb = size(M,1)/NDof;
ii = ((1:NDof)-1)*Nb;
for i = 1:Nb
    D = M(ii + i,ii + i,:);
    Di = myinv(D);
    J(ii + i,ii + i,:) = Di;
end