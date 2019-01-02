function q = rotor_frf(S,R,F,O,w)
NDof = size(S,1)/2;

assert(length(O) == size(S,3));
assert(length(O) == size(R,3));

F = [F; 0*F];

if nargin > 5 %asynchronous
    q = zeros(2*NDof,length(w),length(O));
    assert(length(w) == size(F,2));
    for i = 1:length(O)
        for k = 1:length(w)
            q(:,k,i) = ((1i*w(k)*S(:,:,i)-R(:,:,i))\F(:,k));
        end
    end
    q = q((NDof+1):end,:,:);
else %synchronous so O == w
    q = zeros(2*NDof,length(O));
    assert(length(O) == size(F,2));
    for i = 1:length(O)
        q(:,i) = ((1i*O(i)*S(:,:,i)-R(:,:,i))\F(:,i));
    end
    q = q(NDof+1:end,:,:);
end
