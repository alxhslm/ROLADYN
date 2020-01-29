function [Rr, Ar] = setupFE(Rotor,Bearing)
% SETUPFE performs model reduction on each rotor, including enforcing rigid
% body constraints, and Craig-Bampton reduction

NDofe = 4;

%Handle the rotor DOF first
Ar = [];
Rr = [];
for i = 1:length(Rotor)
    Rcon = zeros(0,Rotor{1}.NDof);
    for j = 1:length(Rotor{i}.Shaft)
        %lock out shaft DOF if shaft is rigid
        SShaft = Rotor{i}.Shaft{j}.S;
        
        if isinf(Rotor{i}.Shaft{j}.Material.E)
            %enforce the displacement of each node to be a
            %rigid body transformations from the previous
            for k = 1:(Rotor{i}.Shaft{j}.Mesh.Nz-1)
                dz = Rotor{i}.Shaft{j}.Element{k}.L;
                Se = Rotor{i}.Shaft{j}.Element{k}.S*SShaft;
                NDofShaft = 4;
                Rcon = [Rcon; axial_offset(dz)*Se(1:NDofShaft,:) - Se(NDofShaft + (1:NDofShaft),:)];
            end
        end
    end
    
    for j = 1:length(Rotor{i}.Disc)
        %lock out disc DOF if shaft is rigid
        SDisc = Rotor{i}.Disc{j}.S;
        SHub  = Rotor{i}.Disc{j}.Hub.S*SDisc;
        SRoot = Rotor{i}.Disc{j}.Root.S*SDisc;
        
        if strcmp(Rotor{i}.Disc{j}.Type,'Flexible')
            if isinf(Rotor{i}.Disc{j}.Material.E)
                %enforce the displacement of each node to be a
                %rigid body transformations from the hub
                for k = 1:Rotor{i}.Disc{j}.Mesh.Nt
                    for l = 1:Rotor{i}.Disc{j}.Mesh.Nr
                        Rcon = [Rcon; Rotor{i}.Disc{j}.Mesh.RHub{k,l}*SHub - Rotor{i}.Disc{j}.Mesh.SNode{k,l}*SDisc];
                    end
                end
            end
            
            %hub
            for k = 1:Rotor{i}.Disc{j}.Mesh.Nt
                SEdge = (Rotor{i}.Disc{j}.Mesh.RHub{k,1}*SHub - Rotor{i}.Disc{j}.Mesh.SNode{k,1}*SDisc);
                
                %axially
                if isinf(Rotor{i}.Disc{j}.Edge.Kzz)
                    Rcon = [Rcon; SEdge(1,:)];
                end
                
                %circumferentially
                if isinf(Rotor{i}.Disc{j}.Edge.Ktt)
                    Rcon = [Rcon; SEdge(2,:)];
                end
                
                %radially
                if isinf(Rotor{i}.Disc{j}.Edge.Krr)
                    Rcon = [Rcon; SEdge(3,:)];
                end
            end
        end
        
        %root
        if isinf(Rotor{i}.Disc{j}.Root.Krr)                
            Rcon = [Rcon; SHub([1 2],:)-SRoot([1 2],:)];
        end
        if isinf(Rotor{i}.Disc{j}.Root.Ktt)                
            Rcon = [Rcon; SHub([3 4],:)-SRoot([3 4],:)];
        end
    end
    
    Kr = Rotor{i}.K;
    for j = 1:length(Rotor{i}.Bearing)
        B = Bearing{Rotor{i}.Bearing{j}.iBearing};
        switch Rotor{i}.Bearing{j}.iNodeBearing
            case 1
                Kb = B.Ri' * B.Kb(1:4,1:4) * B.Ri;
            case 2
                Kb = B.Ro' * B.Kb(5:8,5:8) * B.Ro;
        end
        Sb = Rotor{i}.SNode{Rotor{i}.Bearing{j}.iNode};
        Kr = Kr + Sb'*Kb*Sb;
    end
    
    %apply any rigid shaft/disc constraints
    if isempty(Rcon)
        A = eye(Rotor{i}.NDof);
    else
        A = null(Rcon,'r');
    end
       
    Kr   = A'*Kr*A;
    Mr   = A'*Rotor{i}.M*A;
    Fgr  = A'*Rotor{i}.Fg;
          
    %decide if we want to do Craig-Bampton model reduction
    if isfield(Rotor{i},'iModes')
        %work out which nodes are fixed (to bearings)
        iFixed = [];
        for j = 1:length(Rotor{i}.Bearing)
            if ~Rotor{i}.Bearing{j}.bLinear
                iFixed = [iFixed;
                         (Rotor{i}.Bearing{j}.iNode-1)*NDofe + Rotor{i}.Bearing{j}.iActive];
            end
        end
        
        %compute modes of rotor subsystem
        [Vm,Vc] = cms_analysis(Mr,Kr,Fgr,A,iFixed);
        
        %choose which modes to retain
        if length(Rotor{i}.iModes)>1
            iRetain = Rotor{i}.iModes;
        else
            iRetain = 1:Rotor{i}.iModes;
        end
        
        bKeep = false(size(Vm,2),1);
        bKeep(iRetain) = true;
            
        iOutofBounds = iRetain>size(Vm,2);
        if any(iOutofBounds)
            warning('User specified retaining modes %s which were out of bounds',mat2str(iRetain(iOutofBounds)));
            iRetain(iOutofBounds) = [];
        end
        
        %the reponse must be in the column-space of the retained modes
        Acb = [Vm(:,bKeep) Vc];
        A = A*Acb;
    end
    
    R = null(A')';
    
    %now store everything
    Rotor{i}.FE.A = A;
    Rotor{i}.FE.R = R;
    
    Ar = [Ar Rotor{i}.S'*A];
    Rr = [Rr; R*Rotor{i}.S];
end

function [Vm,Vc] = cms_analysis(M,K,F,A,iFixed)

%perform modal analysis with fixed constraint nodes
Bfree = null(A(iFixed,:));
Bfixed = A(iFixed,:)\eye(length(iFixed));
if isempty(Bfree)
    %no internal modes
    Vm = zeros(size(A,2),0);
    Vc = eye(size(A,2));
    return;
end

%isolate free nodes
Kii = Bfree'*K*Bfree;
Mii = Bfree'*M*Bfree;
Kib = Bfree'*K*Bfixed;
Fi = Bfree'*F;

if ~isempty(iFixed)
    %find constraint modes
    xf = -Kii\Kib;
    Vc = Bfree*xf + Bfixed;
else
    Vc = zeros(size(Bfree,1),0);
end

% %deflection due to self-weight
% Vc(:,end+1) = Bfree*(Kii\Fi);

if rank(Mii) < size(Mii,2)
    B = null(Mii);
    Aslave = null(B'*Kii);
        
    Kii  = Aslave'*Kii*Aslave;
    Mii  = Aslave'*Mii*Aslave;

    Bfree = Bfree*Aslave;
    
    %enforce symmetry
    Mii = 0.5*(Mii + Mii.');
    Kii = 0.5*(Kii + Kii.');
end

[Vfree,D] = eig(Kii,Mii,'vector');
[D,ii] = sort(D);
Vfree = Vfree(:,ii);
Vm = Bfree*Vfree;