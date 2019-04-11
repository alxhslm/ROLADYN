function REB = setupRaceModel(REB)
r = linspace(0,1E-5,100);
r = [r(1:end-1),linspace(1E-5,1E-4,20)];

NDof = REB.Options.bRaceCompliancei + REB.Options.bRaceComplianceo;
if NDof == 0
    REB.Contact.Race.K = REB.Contact.K;
    REB.Contact.Race.n = REB.Contact.n;
    REB.Contact.Race.c = 0;
    return
end

options = optimoptions('fsolve','Display','off','TolX',1E-12,'TolFun',1E-6);
w0 = zeros(NDof,1);
for i = 1:length(r)
    [w(:,i),~,status] = fsolve(@(x)ball_equib(x,r(i),REB),w0,options);
    w0 = w(:,i);
end
[~,Q] = ball_equib(w,r,REB);

if 0%isfield(REB.Contact,'Race')
    p0 = [REB.Contact.Race.K REB.Contact.Race.n REB.Contact.Race.c REB.Contact.e]';
else
    p0 = [REB.Contact.K REB.Contact.n 0 1]';
end
pLB = [0.25*[REB.Contact.K REB.Contact.n] -100 -10]';
pUB = [4*[REB.Contact.K REB.Contact.n]  100 10]';
pscale = [REB.Contact.K REB.Contact.n 10 0.1]';
p = modelfit(r,Q,@poly_hertz,p0./pscale,pLB./pscale,pUB./pscale,pscale);
p = p.*pscale;

REB.Contact.Race.K = p(1);
REB.Contact.Race.n = p(2);
REB.Contact.Race.c = p(3);
REB.Contact.Race.e = p(4);

function Q = poly_hertz(p,r,pscale)
p = p.*pscale;
Q = p(1)*r.^p(2) .* (1 + p(3)*r.^p(4));

function [f,Q_hertz] = ball_equib(x,r,REB)
wi = 0;
wo = 0;
count = 1;
if REB.Options.bRaceCompliancei
    wi = x(count,:);
    count = count + 1;
end
if REB.Options.bRaceComplianceo
    wo = x(count,:);
    count = count + 1;
end
Q_hertz = hertz_contactlaw(REB.Contact.K,REB.Contact.n,r-wo-wi,1E-8);
f = [];

if REB.Options.bRaceCompliancei
    f = [f; Q_hertz - race_compliance(REB.Race.Inner,wi)];
end
if REB.Options.bRaceComplianceo
    f = [f; Q_hertz - race_compliance(REB.Race.Outer,wo)];
end