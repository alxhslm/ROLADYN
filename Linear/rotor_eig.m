function [V,d,W,S,R,A] = rotor_eig(FE,O)
%ROTOT_EIG 

NDof = 2*size(FE.A,2);
NSweep = length(O);

[R,S,A] = rotor_ss(0,FE);
S = repmat(S,1,1,NSweep);
R = repmat(R,1,1,NSweep);
A = repmat(A,1,1,NSweep);

for i = 1:NSweep
    [R(:,:,i),S(:,:,i),A(:,:,i)] = rotor_ss(O(i),FE);
    A2(:,:,i) = blkdiag(A(:,:,i),A(:,:,i));
end

At = mtransposex(A2);
S2 = mtimesx(At,mtimesx(S,A2));
R2 = mtimesx(At,mtimesx(R,A2));

[V,d,W] = eigenshuffle(R2,S2,O);

for i = 1:NSweep
    V2(:,:,i) = A2(:,:,i)*V(:,:,i);
    W2(:,:,i) = A2(:,:,i)*W(:,:,i);
end
V = V2;
W = W2;