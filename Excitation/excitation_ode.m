function [U,Udot,Uddot] = excitation_ode(P,O,w,int_O_dt,int_w_dt,dO_dt,dw_dt)

if nargin < 6
    dO_dt = 0*O;
end
if nargin < 7
    dw_dt = 0*w;
end

U     = zeros(P.Model.Excite.NExcite,length(O));
Udot  = zeros(P.Model.Excite.NExcite,length(O));
Uddot = zeros(P.Model.Excite.NExcite,length(O));

for i = 1:length(P.Excite)
     if strcmpi(P.Excite{i}.Mode, 'Sync')
         phase      = int_O_dt;
         omega      = O;
         domega_dt  = dO_dt;
     elseif strcmpi(P.Excite{i}.Mode, 'Async')
         phase      = int_w_dt;
         omega      = w;
         domega_dt  = dw_dt;
     else
         error('Unknown excitation frequency type');
     end
     
     u     =  P.Excite{i}.u.*exp(1i.*phase);
     udot  =  u*1i.*omega;
     uddot = -u.*omega.^2 + u*1i.*domega_dt;
     
     U     = U   + P.Excite{i}.S'*u;
     Udot  = Udot  + P.Excite{i}.S'*udot;
     Uddot = Uddot + P.Excite{i}.S'*uddot;
end

%remove any singular dimensions
U     = real(U);
Udot  = real(Udot);
Uddot = real(Uddot);
