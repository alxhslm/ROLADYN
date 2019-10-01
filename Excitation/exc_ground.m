function [U,Udot,Uddot] = exc_ground(E,P,w,ph,dw)

U   = zeros(0,numel(w));
Udot  = zeros(0,numel(w));
Uddot = zeros(0,numel(w));

if nargin  < 5
    dw = 0*w;
    if nargin < 4
        ph = 0*w;
    end
end

wons = 0*w + 1;


ph = repmat(ph,4,1);
w  = repmat(w,4,1);
dw = repmat(dw,4,1);

for j = 1:2
    ph  =  1i.*(ph + E.aGndPh{j}*wons);
    u   = (E.uGnd{j}*wons).*exp(ph);
    udot  = (1i.*w).*u;
    uddot = (-w.^2+1i.*dw).*u;

    U     = [U;     u];
    Udot  = [Udot;  udot];
    Uddot = [Uddot; uddot];
end
