function REB = setupDynamicModel(REB)

kc = REB.Dynamics.mb * REB.Geometry.dm/2;
r = unique([linspace(0,15E-6,40) linspace(15E-6,25E-6,10) linspace(25E-6,50E-6,2)]);
Fc = unique([linspace(0,100,20) linspace(100,1000,10)]);
dc = (Fc/REB.Contact.Outer.K).^(1/REB.Contact.Outer.n);
O = sqrt(Fc/kc);
options = optimoptions('fsolve','Display','off','TolX',1E-8,'TolFun',1E-6);
w0 = 0;
for j = 1:length(O)           
    for i = 1:length(r)
        [w(j,i),~,status] = fsolve(@(x)ball_equib(x,r(i),O(j),REB),w0,options);
        w0 = w(j,i);
    end
    w0 = w(j,1);
end

[qr,qO] = meshgrid(r,O);

Qi = hertz_contactlaw(REB.Contact.Inner.K,REB.Contact.Inner.n,qr-w,1E-8);

%fit a model for each speed
K = REB.Contact.K;
n = REB.Contact.n;
% pLB = [1 1 -1 -1 -1];% pLB = p0-Inf;
% pUB = [1 1 -1 1  0];% pUB = p0+Inf;
% p0  = [1 1 -1 0 -0.3];
pLB = [0.25 0.25 -1 ];% pLB = p0-Inf;
pUB = [10    10  -1];% pUB = p0+Inf;
p0  = [1      1  -1 ];

for j = length(O):-1:1
%     pscale = [K n dc(j) 1 1];
    pscale = [K n dc(j)];
    [p(j,:),Qi_fit(j,:)] = modelfit(r,Qi(j,:),@dyn_contact,p0',pLB',pUB',pscale');
    p0 = p(j,:);
    p(j,:) = p(j,:).*pscale;
end

% figure
% surf(r,O,Qi);
% hold on
% surf(r,O,Qi_fit);

REB.Contact.Dynamic = struct('O',O(:),'p',p);
REB.Contact.Dynamic.p  = p;

function Q = dyn_contact(p,r,pscale)
if nargin > 2
p = p.*pscale;
end
r = r + p(3,:);
Q = maxSmooth(p(1,:).*sgn_power(r,p(2,:)),0,1E-8);
% Q = maxSmooth(p(1,:).*sgn_power(r,p(2,:)).*(1+p(4,:).*(abs(r)+eps).^p(5,:)),0,1E-8);

function [f,Q_i,Q_o] = ball_equib(x,r,Ocage,REB)
Fc = (REB.Dynamics.mb * Ocage.^2 * REB.Geometry.dm/2);

Q_i = hertz_contactlaw(REB.Contact.Inner.K,REB.Contact.Inner.n,r-x,1E-8);
Q_o = hertz_contactlaw(REB.Contact.Outer.K,REB.Contact.Outer.n,x,1E-8);

f = Q_i + Fc - Q_o;