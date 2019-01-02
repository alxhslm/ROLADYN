function [FMExcite,uExcite,duExcite,dduExcite] = excitation_ode(P,O,w,int_O_dt,int_w_dt,dO_dt,dw_dt)

NDofTot   = P.FE.NDofFull;
NBearings = length(P.Bearing);
uExcite   = zeros(4*2*NBearings,size(O,1),size(O,2));
duExcite  = zeros(4*2*NBearings,size(O,1),size(O,2));
dduExcite = zeros(4*2*NBearings,size(O,1),size(O,2));
FMExcite  = zeros(NDofTot,size(O,1),size(O,2));

for i = 1:length(P.Excite)
     if strcmpi(P.Excite(i).Mode, 'Sync')
         phase      = int_O_dt;
         omega      = O;
         domega_dt  = dO_dt;
     elseif strcmpi(P.Excite(i).Mode, 'Async')
         phase      = int_w_dt;
         omega      = w;
         domega_dt  = dw_dt;
     else
         error('Unknown excitation frequency type');
     end

     [F,u,du,ddu] = feval(['exc_' P.Excite(i).Name],P.Excite(i).Params,P,omega,phase,domega_dt);
     FMExcite  = FMExcite  + F;
     uExcite   = uExcite   + u;
     duExcite  = duExcite  + du;
     dduExcite = dduExcite + ddu;
end

%remove any singular dimensions
FMExcite  = squeeze(FMExcite);
uExcite   = squeeze(uExcite);
duExcite  = squeeze(duExcite);
dduExcite = squeeze(dduExcite);

FMExcite  = real(FMExcite);
uExcite   = real(uExcite);
duExcite  = real(duExcite);
dduExcite = real(dduExcite);
