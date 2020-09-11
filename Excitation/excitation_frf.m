function [U,Udot,Uddot] = excitation_frf(P,O,w)
if nargin < 3
    w = 0;
end

U = zeros(P.Model.Excite.NExcite,length(O),length(w));
Udot = zeros(P.Model.Excite.NExcite,length(O),length(w));
Uddot = zeros(P.Model.Excite.NExcite,length(O),length(w));

for i = 1:length(P.Excite)
     if strcmpi(P.Excite{i}.Mode, 'Sync')
         omega = O;
     elseif strcmpi(P.Excite{i}.Mode, 'Async')
         omega = w;
     else
         error('Unknown excitation frequency type');
     end
     
     u = P.Excite{i}.u;
     udot = 1i*P.Excite{i}.u*omega;
     uddot = -P.Excite{i}.u*omega.^2;
     
     if strcmpi(P.Excite{i}.Mode, 'Sync')
         u   = repmat(u,1,1,length(w));
         udot  = repmat(udot,1,1,length(w));
         uddot = repmat(uddot,1,1,length(w));
     elseif strcmpi(P.Excite{i}.Mode, 'Async')
         u   = repmat(permute(u,  [1 3 2]),1,length(O),1);
         udot  = repmat(permute(udot, [1 3 2]),1,length(O),1);
         uddot = repmat(permute(uddot,[1 3 2]),1,length(O),1);
     end
     
     U     = U + mtimesx(P.Excite{i}.U',u);
     Udot  = Udot + mtimesx(P.Excite{i}.U',udot);
     Uddot = Uddot + mtimesx(P.Excite{i}.U',uddot);
end