function U = excitation_hbm(P)
U = cell(1,2);

for i = 1:length(P.Excite)
     if strcmpi(P.Excite{i}.Mode, 'Sync')
         iFreq = 1;
     elseif strcmpi(P.Excite{i}.Mode, 'Async')
         iFreq = 2;
     else
         error('Unknown excitation frequency type');
     end
        
     U{iFreq} = [U{iFreq}; P.Excite{i}.u(:)];
end

if isempty(U{2})
    U{2} = 0*U{1};
end