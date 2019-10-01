function [U,Udot,Uddot] = exc_unbalance(E,P,wShaft,aShaft,dwShaft)

if nargin  < 5
    dwShaft = 0.*wShaft;
    if nargin < 4
        aShaft = 0.*wShaft;
    end
end

iRtr = E.iRotor;

a  = aShaft   .* P.Rotor{iRtr}.Speed;
w  = wShaft   .* P.Rotor{iRtr}.Speed;
dw = dwShaft  .* P.Rotor{iRtr}.Speed;

a   = repmat(a,4,1);
w   = repmat(w,4,1);
dw  = repmat(dw,4,1);

wons = a(1,:)*0 + 1;

ph = [E.aUnbalancePh;E.aSkewPh]*wons;
mag = [E.eUnbalance; E.aSkew];

U = mag.*exp(1i.*(a + ph));
Udot = U .* (1i*w);
Uddot = Udot.* (1i*w)+ U .* (1i*dw);