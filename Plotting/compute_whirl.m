function [kappa,ampl] = compute_whirl(P,modes,iNode,options)
sz = size(modes);
if length(sz) < 3
    sz(3) = 1;
end
NSweep = sz(3);
NModes = sz(2);

NNodes = 0;
iNodesIn = [];
iNodesOut = {};
for i = 1:length(P.Rotor)
    NNodes(i) = length(P.Rotor{i}.Nodes);
    iNodesIn = [iNodesIn P.Rotor{i}.iGlobal(1:(4*NNodes(i)))];
    iNodesOut{i} = 1:NNodes(i);
end

NNodesTot = sum(NNodes);

kappa = zeros([NNodes,sz(2:end)]);
ampl =  zeros([NNodes,sz(2:end)]);

modes = mtimesx(P.Model.A(iNodesIn,:),modes);

%we need to scale the modes such that the outer rotor is (arbitrarily) 
%aligned with the real axis which will then allow comparions across all
%speeds
if nargin > 3 && options.bRot
    u = modes(3:4:end,:,:); 
    v = modes(4:4:end,:,:);
else
    u = modes(1:4:end,:,:); 
    v = modes(2:4:end,:,:);
end

%scale by u(iNode) so this will be entirely real and have amplitude 1
scale = repmat(u(iNode,:,:)./(abs(u(iNode,:,:))+eps),NNodesTot,1,1);
scale(scale == 0) = 1;
u = u ./ (scale + eps);
v = v ./ (scale + eps);
u(iNode,:,:) = real(u(iNode,:,:));

if nargin > 3 && options.bDetectZC
    for j = 1:NModes
        ii = 1;
        while true
            del = permute(u(:,j,ii:end),[1 3 2]);
            iZC = any(del(:,1:end-1).*del(:,2:end)< 0 & abs(del(:,1:end-1) - del(:,2:end))> options.tol);
            ii = find(iZC,1)+ii;
            if isempty(ii)
                break;
            end
            u(:,j,ii:end) = -u(:,j,ii:end);
            v(:,j,ii:end) = -v(:,j,ii:end);
        end
    end
end

Vlast = repmat(eye(2),1,1,NNodesTot);
for j = 1:NModes
    for k = 1:NSweep
        for i = 1:NNodesTot
            T = [real(u(i,j,k)) -imag(u(i,j,k));
                 real(v(i,j,k)) -imag(v(i,j,k))];
            H = T*T';
            [V,L] = eig(H,'vector');
            if k > 1
                V = V.*repmat(sign(sum(V .* Vlast(:,:,i),1)),2,1);
            end
            
            phase = wrapToPi(angle(u(i,j,k))-angle(v(i,j,k)));
            %eccentricity parameter, signed by direction of whirl
            kappa(i,j,k) = sqrt(L(1)/(L(2)+eps))*sign(phase+eps);
            
            %there is an issue here that V can flip sign
            %vector of the semi-major axis
            ampl(i,j,k) = sqrt(L(1))*[1 1i]*V(:,1);
            
            Vlast(:,:,i) = V;
        end
        %now sign the amplitude base on the relative phase of u
        ampl(:,j,k) = ampl(:,j,k) .* (1-2*(real(u(:,j,k)) < 0));
    end
end

% %now rescale by the (signed) amplitudes so that we keep smooth transitions
% scale = repmat(ampl(iNode,:,:)./abs(ampl(iNode,:,:)),NNodes,1,1);
% ampl = ampl./(scale + eps);
ampl = real(u);
% scale = repmat(ampl(iNode,:,:),NNodes,1,1);
% ampl = ampl./scale;

%correct sign of kappa depending on rotor spin direction
for i = 1:length(P.Rotor)
    kappa(iNodesOut{i},:,:,:,:) = kappa(iNodesOut{i},:,:,:,:) * sign(P.Rotor{i}.Speed);
end