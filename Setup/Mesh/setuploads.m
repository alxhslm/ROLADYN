function P = setuploads(P)
NDofTot = P.Mesh.NDof;

Fc = -(P.Mesh.R) \ (P.Mesh.K * P.Mesh.x0  - P.Mesh.Fg);

%% Begin with the rotors
NRotor = length(P.Rotor);
F0 = zeros(NDofTot,1);
for i = 1:NRotor
    P.Rotor{i}.F0 = zeros(P.Rotor{i}.NDof,1);
    %deal with the shafts first
    for j = 1:length(P.Rotor{i}.Shaft)
        SShaft = P.Rotor{i}.Shaft{j}.S;
        
        for k = 1:(P.Rotor{i}.Shaft{j}.Mesh.Nz-1)
            %rotor->shaft element mapping matrix
            Se = P.Rotor{i}.Shaft{j}.Element{k}.S*SShaft;
            Re = P.Rotor{i}.Shaft{j}.Element{k}.R;
            P.Rotor{i}.F0 = P.Rotor{i}.F0 + Se'*Re'*P.Rotor{i}.Shaft{j}.Element{k}.F0;
        end
    end
    
    %now deal with the discs
    for j = 1:length(P.Rotor{i}.Disc)
        %create rotor->disc mapping matrix
        SDisc = P.Rotor{i}.Disc{j}.S;
        
        SHub = P.Rotor{i}.Disc{j}.Hub.S*SDisc;
        SRoot = P.Rotor{i}.Disc{j}.Root.S*SDisc;
        
        if strcmp(P.Rotor{i}.Disc{j}.Type,'Flexible')
            %disc compliance
            for k = 1:P.Rotor{i}.Disc{j}.Mesh.Nt
                for l = 1:(P.Rotor{i}.Disc{j}.Mesh.Nr-1)
                    P.Rotor{i}.F0 = P.Rotor{i}.F0 + Se'*Re'*P.Rotor{i}.Disc{j}.Element{k,l}.F0;
                end
            end
        end
        
        %lump inertia at hub
        P.Rotor{i}.F0 = P.Rotor{i}.F0 + (SHub-SRoot)'*P.Rotor{i}.Disc{j}.Root.F0;
    end
    F0 = F0 + P.Rotor{i}.S'*P.Rotor{i}.F0;
end

P.Mesh.Rotor.F0 = F0 + P.Mesh.Rotor.R*Fc;

%% And the stator
NStator = length(P.Stator);
F0 = zeros(NDofTot,1);
for i = 1:NStator
    F0 = F0 + P.Stator{i}.S'*P.Stator{i}.F0;
end
P.Mesh.Stator.F0 = F0 + P.Mesh.Stator.R*Fc;

%% Move onto the bearings
NBearings = length(P.Bearing);
F0 = zeros(NDofTot,1);
F0b = zeros(2*4*NBearings,1);
for i = 1:NBearings
    F0 = F0 + P.Bearing{i}.S'*P.Bearing{i}.R'*P.Bearing{i}.Fb;
    F0b = F0b + [P.Bearing{i}.Ui;P.Bearing{i}.Uo]'*P.Bearing{i}.R'*P.Bearing{i}.Fb;
end
P.Mesh.Bearing.F0 = F0  + P.Mesh.Bearing.R*Fc;
P.Mesh.Bearing.Fb = F0b + P.Mesh.Bearing.Rb*Fc;

%% Combined
P.Mesh.F0 = P.Mesh.Rotor.F0 + P.Mesh.Bearing.F0 + P.Mesh.Stator.F0;