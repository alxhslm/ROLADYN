function [W,Ferr,V,S] = REB_harris_analytical(B,States)

q    = States.qi - States.qo;
xInt = States.xInt;

Oi = States.Oi;
Oo = States.Oo;
Ai = States.Ai;
Ao = States.Ao;

%create some vectors of the right size
wons = (B.psi*0+1);
x0 = (q(1,:)*0 + 1);

%if we're not including VC effects, set Ocage to 0 to "fix" the cage
if ~B.bVC
    Acage = 0;
else
    Acage = B.rCagei * Ai + B.rCageo * Ao;
end

z = B.zi*sign(B.z/(B.z(1)+eps));
psi = B.psi + Acage;
Z = z*x0;
PSI = psi*x0;

dz  = wons*q(3,:) + B.Ri*(sin(PSI).*(wons*q(4,:)) - cos(PSI).*(wons*q(5,:))) - B.cz;
dr  = cos(PSI).*(wons*q(1,:)) + sin(PSI).*(wons*q(2,:)) - Z.*(sin(PSI).*(wons*q(4,:)) - cos(PSI).*(wons*q(5,:))) - B.cr;

Az = B.A0*sin(B.alpha)*x0 + dz;
Ar = B.A0*cos(B.alpha)*x0 + dr;
A = sqrt(Az.^2 + Ar.^2);

%now find the ball forces
if B.bCentrifugal || B.bGyro
    vz = xInt(1:B.Z,:);
    vr = xInt(B.Z+(1:B.Z),:);
    Xz = (B.ro-B.D/2)*sin(B.alpha)*x0 + vz;
    Xr = (B.ro-B.D/2)*cos(B.alpha)*x0 + vr;
    [Ai,Ao,alpha_i,alpha_o] = race_geometry(Xz,Xr,Az,Ar);
    dbi = Ai - (B.ri-B.D/2);
    dbo = Ao - (B.ro-B.D/2);
    if nargout > 8
        [Qi,Ki] = hertz_contactlaw(B.Inner.K,B.n,dbi);
        [Qo,Ko] = hertz_contactlaw(B.Outer.K,B.n,dbo);
        [Jb,Juu,Juv,Jvu,Jvv] = feval(['REB_harris_' B.Control],Ki,Ko,Ai,Ao,alpha_i,alpha_o,Qi,Qo,wons*Oi,wons*Oo,B);
    else
        Qi = hertz_contactlaw(B.Inner.K,B.n,dbi);
        Qo = hertz_contactlaw(B.Outer.K,B.n,dbo);
    end
    if 0
        [Fc,Fi,Fo,Mg] = dynamic_ball_loads(B,alpha_i,alpha_o,wons*Oi,wons*Oo);
    else
        Fc = 0.5*B.dm*B.mBall*(wons*(B.rCagei*Oi + (1-B.rCageo)*Oo)).^2;
        Fi = 0*Qi;
        Fo = 0*Qo;
        Mg = 0*wons;
    end
    Gu = cat(3,Qi.*sin(alpha_i) - Fi.*cos(alpha_i), Qi.*cos(alpha_i) + Fi.*sin(alpha_i));
    Gv = cat(3,Qi.*sin(alpha_i) - Fi.*cos(alpha_i) - Qo.*sin(alpha_o) + Fo.*cos(alpha_o), ...
               Qi.*cos(alpha_i) + Fi.*sin(alpha_i) - Qo.*cos(alpha_o) - Fo.*sin(alpha_o) + Fc);
else
    alpha = atan2(Az,Ar);
    db = A - B.A0;
    if nargout > 8
        [Q,K] = hertz_contactlaw(B.K,B.n,db);
        Jb = REB_harris_simple(K,A,alpha,Q);
        Juu = Jb;
        Juv = zeros([size(K) 2 0]); 
        Jvu = zeros([size(K) 0 2]);
        Jvv = zeros([size(K) 0 0]);
    else
        Q = hertz_contactlaw(B.K,B.n,db);
    end
    Gu = cat(3,Q.*sin(alpha), Q.*cos(alpha));
    Gv = zeros([size(Q),0]);
    
    lambda = (B.Outer.K / B.Inner.K)^(1/B.n);
    Xz = (Az./A) .* ((B.ro-B.D/2) + max(A - B.A0,0) /(1 + lambda));
    Xr = (Ar./A) .* ((B.ro-B.D/2) + max(A - B.A0,0) /(1 + lambda));
    Qi = Q;
    Qo = Q;
    alpha_i = alpha;
    alpha_o = alpha;
end

Gu  = permute(Gu,[3,4,2,1]);
Gv  = permute(Gv,[3,4,2,1]);
W = zeros(6,1,size(q,2));
Ferr = zeros(B.NDof*B.Z,1,size(q,2));

if nargout < 9
    return;
end

Jb = permute(Jb,[3,4,2,1]);
Kb = zeros(6,6,size(q,2));

if nargout > 3
    Juu = permute(Juu,[3,4,2,1]);
    Juv = permute(Juv,[3,4,2,1]);
    Jvu = permute(Jvu,[3,4,2,1]);
    Jvv = permute(Jvv,[3,4,2,1]);

    Kuu = zeros(6,6,size(q,2));
    Kuv = zeros(6,B.NDof*B.Z,size(q,2));
    Kvu = zeros(B.NDof*B.Z,6,size(q,2));
    Kvv = zeros(B.NDof*B.Z,B.NDof*B.Z,size(q,2));
%     else %want the stiffness
%         Kuu = zeros(6,6,size(q,2));
%         for j = 1:B.Z
%             Juu(:,:,:,j) = Juu(:,:,:,j) - mmult(Juv(:,:,:,j),mmult(minv(Jvv(:,:,:,j)),Jvu(:,:,:,j)));
%         end
end


Z   = permute(Z,[3 4 2 1]);
PSI = permute(PSI,[3 4 2 1]);
c = cos(PSI);
s = sin(PSI);
zero = 0*PSI;
one = zero+1;
R = [  zero       c
       zero       s
       one       zero
    B.Ri*s      -Z.*s
   -B.Ri*c       Z.*c
      zero       zero];
Rt = permute(R,[2 1 3 4]);
   
for j = 1:B.Z
    W = W + mtimesx(R(:,:,:,j),Gu(:,:,:,j));
    for k = 1:B.NDof
        Ferr((k-1)*B.Z+j,:) = Gv(k,:,:,j);
    end
    if nargout > 8
        Kb = Kb + mtimesx(R(:,:,:,j),mtimesx(Jb(:,:,:,j),Rt(:,:,:,j)));
        if nargout > 9
            Kuu = Kuu + mtimesx(R(:,:,:,j),mtimesx(Juu(:,:,:,j),Rt(:,:,:,j)));
            for k = 1:B.NDof
                Kuv(:,(k-1)*B.Z+j,:) = Kuv(:,(k-1)*B.Z+j,:) + mtimesx(R(:,:,:,j),Juv(:,k,:,j));
                Kvu((k-1)*B.Z+j,:,:) = Kvu((k-1)*B.Z+j,:,:) + mtimesx(Jvu(k,:,:,j),Rt(:,:,:,j));
                for l = 1:B.NDof
                    Kvv((k-1)*B.Z+j,(l-1)*B.Z+j,:) = Jvv(k,l,:,j);
                end
            end
        end
    end
end

W = permute(W,[1 3 4 2]);
Ferr = permute(Ferr,[1 3 4 2]);

%contact loads
V.Qi = Qi; V.Qo = Qo; 
V.Fi = Fi; V.Fo = Fo;

%geometry
V.alpha_i = alpha_i; V.alpha_o = alpha_o;
V.dbi = dbi; V.dbo = dbo;
V.Ai = Ai;  V.Ao = Ao;
V.Xr = Xr;  
V.Xz = Xz;

%dynamic loads
V.Fc = Fc;   
V.Mg = Mg;
 
%stiffnesses
S.K = Kb;
S.Jqq = Kuu;
S.Jqx = Kuv;
S.Jxq = Kvu;
S.Jxx = Kvv;

function C = mmult(A,B)
if size(A,3) < 1
    C = A*B;
    return;
end
if isempty(A) || isempty(B)
    C = 0;
    return;
end
C = zeros(size(A,1),size(B,2),size(A,3));

for i = 1:size(A,1)
    for j = 1:size(B,2)
        C(i,j,:) = sum(A(i,:,:).*permute(B(:,j,:),[2 1 3]),2);
    end
end

function A = minv(A)
if isempty(A)
    A = 0;
    return;
end
a = A(1,1,:);
b = A(1,2,:);
c = A(2,1,:);
d = A(2,2,:);
det = a.*d - b.*c;
a = a ./ det;
b = b ./ det;
c = c ./ det;
d = d ./ det;
A = [d -b;
    -c  a];