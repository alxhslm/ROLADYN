function [Forces,Channels,Damping] = SFD_model(Params, States)
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

model = Params.fun;
NPts = size(States.qi,2);

States = default_inputs(States,R);
States = default_speeds(States,NPts);

if nargout < 3
    [Forces,Channels] = feval(model,Params,States);
else
    [Forces,Channels,Damping] = feval(model,Params,States);

    F0 = Forces.F;
    qsgn = sign((States.qi - States.qo)+eps);

    if ~isfield(Damping,'K')
        h = 1E-8;
        KBearing = zeros(6,6,NPts);
        q = States.qi;
        for i = 1:6
            hSigned = -h*qsgn(i,:);
            States.qi(i,:) = States.qi(i,:) + h*hSigned;
            Forces = feval(model,Params,States);
            KBearing(:,i,:) = permute((Forces.F-F0)./(ones(6,1)*hSigned),[1 3 2]);
            States.qi = q; 
        end
        Damping.K = KBearing; 
    end

    if ~isfield(Damping,'C')
        h = 1E-8;
        CBearing = zeros(6,6,NPts);
        qdot = States.qidot;
        for i = 1:6
            hSigned = -h*qsgn(i,:);
            States.qidot(i,:) = States.qidot(i,:) + hSigned;
            Forces = feval(model,Params,States);
            CBearing(:,i,:) = permute((Forces.F-F0)./(ones(6,1)*hSigned),[1 3 2]);
            States.qidot = qdot; 
        end
        Damping.C = CBearing;
    end

    if ~isfield(Damping,'M')
        h = 1E-8;
        MBearing = zeros(6,6,NPts);
        qddot = States.qiddot;
        for i = 1:6
            hSigned = -h*qsgn(i,:);
            States.qiddot(i,:) = States.qiddot(i,:) + hSigned;
            Forces = feval(model,Params,States);          
            MBearing(:,i,:) = permute((Forces.F-F0)./(ones(6,1)*hSigned),[1 3 2]);
            States.qiddot = qddot; 
        end
        Damping.M = MBearing;
    end
end

Forces.F  = kron([1; -1],R'*Forces.F);

Forces.q = R'*(States.qi - States.qo);

if nargout > 2
    %Combine stiffness terms
    Damping.K = kron([1 -1; -1 1],mtimesx(R',mtimesx(Damping.K,R)));
    Damping.C = kron([1 -1; -1 1],mtimesx(R',mtimesx(Damping.C,R)));
    Damping.M = kron([1 -1; -1 1],mtimesx(R',mtimesx(Damping.M,R)));
    
    Damping = default_damping(Damping);
end

function Damping = default_damping(Damping)
field = 'KCM';
for i = 1:3
    if ~isfield(Damping,[field(i) 'qq'])
        K = Damping.(field(i));
        Damping.([field(i) 'qq']) = K;
        Damping.([field(i) 'qx']) = zeros(size(K,1),0,size(K,3));
        Damping.([field(i) 'xq']) = zeros(0,size(K,2),size(K,3));
        Damping.([field(i) 'xx']) = zeros(0,0,size(K,3));
    end
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