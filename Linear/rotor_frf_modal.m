function q = rotor_frf_modal(V,d,W,F,O,w)
NModes = size(d,1);
NDof = size(V,1)/2;

assert(length(O) == size(V,3));
assert(length(O) == size(d,2));
assert(length(O) == size(W,3));

F = [F; 0*F];

if nargin > 5 %asynchronous
    q = zeros(2*NDof,length(w),length(O));
    assert(length(w) == size(F,2));
    for i = 1:length(O)
        for k = 1:NModes
            for j = 1:length(w)
                q(:,j,i) = q(:,j,i) + ((V(:,k,i)*W(:,k,i)')/(1i*w(j) - d(k,i))) * F(:,j);
            end
        end
    end
else %synchronous so O == w
    q = zeros(2*NDof,length(O));
    assert(length(O) == size(F,2));
    for i = 1:length(O)
        for k = 1:NModes
            q(:,i) = q(:,i) +  ((V(:,k,i)*W(:,k,i)')/(1i*O(i) - d(k,i))) * F(:,i);
        end   
    end
end
q = q(NDof+1:end,:,:);