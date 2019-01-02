function [omega,zeta,modes] = eig2modes(V,d,iSort)

NModes = size(d,1)/2;
NSweep = size(d,2);
NDof   = size(V,1)/2;

omega = zeros(NModes,NSweep);
zeta  = zeros(NModes,NSweep);
modes = zeros(NDof,NModes,NSweep);

%group by the complex parts into complex conjugate pairs
[~,ii] = sort(abs(imag(d(:,1))));
d = d(ii,:);
V = V(:,ii,:);

%now ensure that the +ve imag part comes first
for j = 1:size(d,2)
    for i = 1:2:size(d,1)
        if abs(imag(d(i,j)))>1E-10 && imag(d(i,j)) < 0
           d(i+[0 1],j) = d(i+[1 0],j);
           V(:,i+[0 1],j) = V(:,i+[1 0],j);
        end
    end
end

for i = 1:NSweep
    e = d(:,i);
    
    %extract real and imaginary parts
    mu = 0.5*(e(1:2:end) + e(2:2:end));
    sig = 0.5*(e(1:2:end) - e(2:2:end));
    
    %convert real + imag -> omega & zeta
    omega(:,i) = sqrt(mu.^2 - sig.^2);
    zeta(:,i) = -mu./omega(:,i);
    
    %trim off to just the displacement DOF
    modes(:,:,i) = V(NDof+1:end,1:2:end,i);
end

omega = real(omega);
zeta = real(zeta);

%sort if necessary
if nargin > 2   
    [~,ii] = sort(omega(:,iSort));
    omega  = omega(ii,:);
    zeta   = zeta(ii,:);
    modes  = modes(:,ii,:);
end