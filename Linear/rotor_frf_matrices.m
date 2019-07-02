function q = rotor_frf_matrices(K,C,M,F,O,w)
NDof = size(K,1)/2;

assert(length(O) == size(M,3));
assert(length(O) == size(C,3));
assert(length(O) == size(K,3));

if nargin > 5 %asynchronous
    q = zeros(NDof,length(w),length(O));
    assert(length(w) == size(F,2));
    for i = 1:length(O)
        for k = 1:length(w)
            q(:,k,i) = ((K(:,:,i) + 1i*w(k)*C(:,:,i)-w(k)^2*M(:,:,i))\F(:,k));
        end
    end
else %synchronous so O == w
    q = zeros(2*NDof,length(O));
    assert(length(O) == size(F,2));
    for i = 1:length(O)
        q(:,i) = ((K(:,:,i) + 1i*O(i)*C(:,:,i)-O(i)^2*M(:,:,i))\F(:,i));
    end
end
