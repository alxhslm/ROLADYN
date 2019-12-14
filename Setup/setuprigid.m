function [Rr,Ar] = setuprigid(Rotor)
NDofe = 4;

%Handle the rotor DOF first
Ar = [];
Rr = [];
for i = 1:length(Rotor)
    A = zeros(Rotor{i}.NDof,NDofe);
    for j = 1:length(Rotor{i}.Nodes)
        %rigid body constraint
        A((j - 1)*4 + (1:4),:) = axial_offset(Rotor{i}.Nodes(j) - Rotor{i}.z);
    end
    
    %rigid body constraints for discs   
    for j = 1:length(Rotor{i}.Disc)
        %lock out disc DOF if shaft is rigid
        SDisc = Rotor{i}.Disc{j}.S;
        SHub  = Rotor{i}.Disc{j}.Hub.S*SDisc;
        SRoot = Rotor{i}.Disc{j}.Root.S*SDisc;
        
        D = Rotor{i}.Disc{j};
        
        %hub fixed to root
        A = A + SHub'*(SRoot*A);
        
        if strcmp(Rotor{i}.Disc{j}.Type,'Flexible')
            %enforce the displacement of each node to be a
            %rigid body transformations from the hub
            for k = 1:Rotor{i}.Disc{j}.Mesh.Nt
                for l = 1:Rotor{i}.Disc{j}.Mesh.Nr
                    A = A + (D.Mesh.SNode{k,l}*SDisc)'*(D.Mesh.RHub{k,l}*SHub*A);
                end
            end
        end
    end
    
     R = null(A')';
     
     Rotor{i}.Rigid.A = A;
     Rotor{i}.Rigid.R = R;
     
     Ar = [Ar Rotor{i}.S'*A];
     Rr = [Rr; R*Rotor{i}.S];
end