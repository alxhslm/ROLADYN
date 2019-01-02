function [Vseq,Dseq,Wseq] = eigenshuffle(varargin)
% eigenshuffle: Consistent sorting for an eigenvalue/vector sequence
% [Vseq,Dseq] = eigenshuffle(Asequence)
%
% Includes munkres.m (by gracious permission from Yi Cao)
% to choose the appropriate permutation. This greatly
% enhances the speed of eigenshuffle over my previous
% release.
%
% http://www.mathworks.com/matlabcentral/fileexchange/20652
%
% Arguments: (input)
%  Asequence - an array of eigenvalue problems. If
%      Asequence is a 3-d numeric array, then each
%      plane of Asequence must contain a square
%      matrix that will be used to call eig.
%
%      Eig will be called on each of these matrices
%      to produce a series of eigenvalues/vectors,
%      one such set for each eigenvalue problem.
%
% Arguments: (Output)
%  Vseq - a 3-d array (pxpxn) of eigenvectors. Each
%      plane of the array will be sorted into a
%      consistent order with the other eigenvalue
%      problems. The ordering chosen will be one
%      that maximizes the energy of the consecutive
%      eigensystems relative to each other.
%
%  Dseq - pxn array of eigen values, sorted in order
%      to be consistent with each other and with the
%      eigenvectors in Vseq.
%
% Example:
%  Efun = @(t) [1 2*t+1 t^2 t^3;2*t+1 2-t t^2 1-t^3; ...
%               t^2 t^2 3-2*t t^2;t^3 1-t^3 t^2 4-3*t];
%
%  Aseq = zeros(4,4,21);
%  for i = 1:21
%    Aseq(:,:,i) = Efun((i-11)/10);
%  end
%  [Vseq,Dseq] = eigenshuffle(Aseq);
%  
% To see that eigenshuffle has done its work correctly,
% look at the eigenvalues in sequence, after the shuffle.
%
% t = (-1:.1:1)';
% [t,Dseq']
% ans =
%        -1     8.4535           5      2.3447     0.20181
%      -0.9     7.8121      4.7687      2.3728     0.44644
%      -0.8     7.2481        4.56      2.3413     0.65054
%      -0.7     6.7524      4.3648      2.2709      0.8118
%      -0.6     6.3156      4.1751      2.1857     0.92364
%      -0.5     5.9283      3.9855      2.1118     0.97445
%      -0.4     5.5816      3.7931      2.0727     0.95254
%      -0.3     5.2676      3.5976      2.0768       0.858
%      -0.2     4.9791      3.3995      2.1156     0.70581
%      -0.1     4.7109         3.2      2.1742     0.51494
%         0     4.4605           3      2.2391     0.30037
%       0.1     4.2302         2.8      2.2971    0.072689
%       0.2     4.0303      2.5997      2.3303    -0.16034
%       0.3     3.8817      2.4047      2.3064    -0.39272
%       0.4     3.8108      2.1464      2.2628    -0.62001
%       0.5     3.8302      1.8986      2.1111    -0.83992
%       0.6     3.9301      1.5937      1.9298     -1.0537
%       0.7     4.0927      1.2308       1.745     -1.2685
%       0.8     4.3042     0.82515      1.5729     -1.5023
%       0.9     4.5572     0.40389      1.4272     -1.7883
%         1     4.8482  -8.0012e-16     1.3273     -2.1755
%
% Here, the columns are the shuffled eigenvalues.
% See that the second eigenvalue goes to zero, but
% the third eigenvalue remains positive. We can plot
% eigenvalues and see that they have crossed, near
% t = 0.35 in Efun.
%
% plot(-1:.1:1,Dseq')
%
% For a better appreciation of what eigenshuffle did,
% compare the result of eig directly on Efun(.3) and
% Efun(.4). Thus:
%
% [V3,D3] = eig(Efun(.3))
% V3 =
%     -0.74139      0.53464     -0.23551       0.3302
%      0.64781       0.4706     -0.16256      0.57659
%    0.0086542     -0.44236     -0.89119      0.10006
%     -0.17496     -0.54498      0.35197      0.74061
%
% D3 =
%     -0.39272            0            0            0
%            0       2.3064            0            0
%            0            0       2.4047            0
%            0            0            0       3.8817
%
% [V4,D4] = eig(Efun(.4))
% V4 =
%     -0.73026      0.19752      0.49743      0.42459
%      0.66202      0.21373      0.35297      0.62567
%     0.013412     -0.95225      0.25513      0.16717
%     -0.16815    -0.092308     -0.75026      0.63271
%
% D4 =
%     -0.62001            0            0            0
%            0       2.1464            0            0
%            0            0       2.2628            0
%            0            0            0       3.8108
%
% With no sort or shuffle applied, look at V3(:,3). See
% that it is really closest to V4(:,2), but with a sign
% flip. Since the signs on the eigenvectors are arbitrary,
% the sign is changed, and the most consistent sequence
% will be chosen. By way of comparison, see how the
% eigenvectors in Vseq have been shuffled, the signs
% swapped appropriately.
%
% Vseq(:,:,14)
% ans =
%       0.3302      0.23551     -0.53464      0.74139
%      0.57659      0.16256      -0.4706     -0.64781
%      0.10006      0.89119      0.44236   -0.0086542
%      0.74061     -0.35197      0.54498      0.17496
%
% Vseq(:,:,15)
% ans =
%      0.42459     -0.19752     -0.49743      0.73026
%      0.62567     -0.21373     -0.35297     -0.66202
%      0.16717      0.95225     -0.25513    -0.013412
%      0.63271     0.092308      0.75026      0.16815
%
% See also: eig
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 3.0
% Release date: 2/18/09

if nargin < 2
    error('Not enough input arguments');
elseif nargin < 3
    Aseq = varargin{1};
    Bseq = eye(size(Aseq,1));
    O = varargin{2};
else
    Aseq = varargin{1};
    Bseq = varargin{2};
    O = varargin{3};
end

% Is Asequence a 3-d array?
Asize = size(Aseq);
Bsize = size(Bseq);
   

if length(Asize)<3
    Asize(3) = 1;
end

if length(Bsize)<3
    Bsize(3) = 1;
end

if (Asize(1)~=Asize(2))
    error('Aseq must be a (pxpxn) array of eigen-problems, each of size pxp')
end
if  (Bsize(1)~=Bsize(2))
	error('Bseq must be a (pxpxn) array of eigen-problems, each of size pxp')
end
if (Asize(1)~=Bsize(1))
	error('Aseq and Bsequences must both be (pxpxn) arrays of the same size p')
end
if (Asize(3)~=Bsize(3))
    error('Aseq and Bsequences must both be (pxpxn) arrays with the same number of eigen-problems n')
end

p = Asize(1);
if length(Asize)<3
  n = 1;
else
  n = Asize(3);
end

% the initial eigenvalues/vectors in nominal order
Vseq = zeros(p,p,n);
Wseq = zeros(p,p,n);
Dseq = zeros(p,n);
for i = 1:n
  [V,D,W] = eig(Aseq(:,:,i),Bseq(:,:,i),'vector');
%   [V,D,W] = eig(Bseq(:,:,i)\Aseq(:,:,i),'vector');
  % initial ordering is purely in decreasing order.
  % If any are complex, the sort is in terms of the
  % real part.
  [~,iSort] = sort(real(D),1,'descend');
  
  [V,W] = mkorthog(Bseq(:,:,i),V,W);
   
  chkorthog(Aseq(:,:,i),Bseq(:,:,i),V,D,W)
    
  Dseq(:,i) = D(iSort);
  Vseq(:,:,i) = V(:,iSort);
  Wseq(:,:,i) = W(:,iSort);
end

% was there only one eigenvalue problem?
if n < 2
  % we can quit now, having sorted the eigenvalues
  % as best as we could.
  return
end

% now, treat each eigenproblem in sequence (after
% the first one.)
for i = 2:n
  % compute distance between systems
  
  V2 = Vseq(:,:,i);
  D2 = Dseq(:,i);
  W2 = Wseq(:,:,i);
  
  %now predict the new W and D by extrapolating from the old ones
  if i > 2
      dVdO = (Vseq(:,:,i-1) - Vseq(:,:,i-2)) / (O(i-1) - O(i-2) + eps);
      dDdO = (Dseq(:,i-1)   - Dseq(:,i-2)  ) / (O(i-1) - O(i-2) + eps);
      dWdO = (Wseq(:,:,i-1) - Wseq(:,:,i-2)) / (O(i-1) - O(i-2) + eps);
  else
      dVdO = zeros(p,p);
      dDdO = zeros(p,1);
      dWdO = zeros(p,p);
  end
  
  V1 = Vseq(:,:,i-1) + dVdO * (O(i) - O(i-1));
  D1 = Dseq(:,i-1)   + dDdO * (O(i) - O(i-1));
  W1 = Wseq(:,:,i-1) + dWdO * (O(i) - O(i-1));
  
  %Penalise errors in both the eigenvectors and -values
  dist_vec = abs(1-W1'*Bseq(:,:,i)*V2);
  dist_val = sqrt(distancematrix(real(D1),real(D2)).^2+ ...
                  distancematrix(imag(D1),imag(D2)).^2);
  
  % Is there a best permutation? use munkres.
  % much faster than my own mintrace, munkres
  % is used by gracious permission from Yi Cao.
  reorder = munkres(dist_vec.*dist_val);
  
  V2 = V2(:,reorder);
  D2 = D2(reorder);
  W2 = W2(:,reorder);

  %if we have any conjugate pair eigenvalues mixed up
  P = find(imag(D2).* imag(D1)<0);
  if ~isempty(P)
      P = P(1:length(P) - mod(length(P),2));
      P2 = reshape(flipud(reshape(P,2,[])),[],1);
      V2(:,P) = V2(:,P2);
      W2(:,P) = W2(:,P2);
      D2(P) = D2(P2);
  end
  
   % also ensure the signs of each eigenvector pair
   % were consistent if possible
   S = squeeze(real(sum(V1.*V2,1))) < 0;
   V2(:,S) = -V2(:,S);
   W2(:,S) = -W2(:,S);
   
   Vseq(:,:,i) = V2;
   Dseq(:,i)   = D2;
   Wseq(:,:,i) = W2;
    
   chkorthog(Aseq(:,:,i),Bseq(:,:,i),Vseq(:,:,i),Dseq(:,i),Wseq(:,:,i));
end

function chkorthog(A,B,V,D,W)
tol = 1E-3;
err = abs(diag(W'*B*V) - 1);
% assert(all(err<tol));

err = abs(diag(W'*A*V) - D)./(abs(D)+eps);
% assert(all(err(~isinf(D) & D > 1E-3)<tol));

function [V,W] = mkorthog(B,V,W)
 %ensure we normalise the modes correctly
scale = sqrt(diag(W' * B * V)).';    
scale(scale == 0) = 1;

%now scale the eigenvectors
p = size(V,1);
V = V./ repmat(scale,p,1);
W = W./ repmat(conj(scale),p,1);

% =================
% end mainline
% =================
% begin subfunctions
% =================

function d = distancematrix(vec1,vec2)
% simple interpoint distance matrix
[vec1,vec2] = ndgrid(vec1,vec2);
d = abs(vec1 - vec2);
