function [U,Udot,Uddot] = excitation_frf(P,O,w)
if nargin < 3
    w = 0;
end

U = zeros(0,length(O),length(w));
Udot = zeros(0,length(O),length(w));
Uddot = zeros(0,length(O),length(w));

for i = 1:length(P.Excite)
     if strcmpi(P.Excite{i}.Mode, 'Sync')
         omega = O;
     elseif strcmpi(P.Excite{i}.Mode, 'Async')
         omega = w;
     else
         error('Unknown excitation frequency type');
     end

     [u,udot,uddot] = feval(['exc_' P.Excite{i}.Name],P.Excite{i},P,omega);
      if strcmpi(P.Excite{i}.Mode, 'Sync')
         u   = repmat(u,1,1,length(w));
         udot  = repmat(udot,1,1,length(w));
         uddot = repmat(uddot,1,1,length(w));
      elseif strcmpi(P.Excite{i}.Mode, 'Async')
         u   = repmat(permute(u,  [1 3 2]),1,length(O),1);
         udot  = repmat(permute(udot, [1 3 2]),1,length(O),1);
         uddot = repmat(permute(uddot,[1 3 2]),1,length(O),1);
      end
      
      U     = [U;     u];
      Udot  = [Udot;  udot];
      Uddot = [Uddot; uddot];
end