function P = setupmatrices(P,part)

info.NDofe = 4; %4 dof per elem
info.gvec = [P.g; 0; 0];
info.IMap = eye(P.Mesh.NDof);
info.NDof = P.Mesh.NDof;
info.NDofBearing = P.Mesh.Bearing.NDof;
info.NDofRotor = P.Mesh.Rotor.NDof;

if ~iscell(part)
    part = {part};
end

for k = 1:length(part)
    switch part{k}
        case 'Rotor'
            [P.Rotor,P.Mesh.Rotor] = setuprotors(P.Rotor,P.Mesh.Rotor,info);
        case 'Bearing'
            [P.Bearing,P.Mesh.Bearing] = setupbearings(P.Bearing,P.Mesh.Bearing,P.Rotor,info);
        case 'Excite'
            [P.Excite,P.Mesh.Excite] = setupexcitation(P.Excite,P.Mesh.Excite,P.Rotor,P.Bearing,info);
    end
end

function [Rotor,Mesh] = setuprotors(Rotor,Mesh,info)

Mesh.K = zeros(info.NDof);
Mesh.C = zeros(info.NDof);
Mesh.M = zeros(info.NDof);
Mesh.G = zeros(info.NDof);
Mesh.F0 = zeros(info.NDof,1);
Mesh.Fg = zeros(info.NDof,1);

for i = 1:length(Rotor)
    %Mesh.Ceate the global -> rotor mapping matrix
    Sr  = info.IMap(Rotor{i}.iGlobal,:);
    Rotor{i}.S = Sr;
        
    Rotor{i}.M = zeros(Rotor{i}.NDof);
    Rotor{i}.G = zeros(Rotor{i}.NDof);
    Rotor{i}.C = zeros(Rotor{i}.NDof);
    Rotor{i}.K = zeros(Rotor{i}.NDof);
    Rotor{i}.Fg = zeros(Rotor{i}.NDof,1);
    Rotor{i}.F0 = zeros(Rotor{i}.NDof,1);
    
    IMapRotor = eye(Rotor{i}.NDof);
    
    %deal with the shafts first
    for j = 1:length(Rotor{i}.Shaft)
        SShaft = IMapRotor(Rotor{i}.Shaft{j}.iLocal,:);
        Rotor{i}.Shaft{j}.S = SShaft;
        
        for k = 1:(Rotor{i}.Shaft{j}.Nz-1)
            %rotor->shaft element mapping matrix
            Se = Rotor{i}.Shaft{j}.Se{k}*SShaft;
                        
            %stiffness/damping
            Rotor{i}.K = Rotor{i}.K + Se'*Rotor{i}.Shaft{j}.R{k}'*Rotor{i}.Shaft{j}.Ke{k}*Rotor{i}.Shaft{j}.R{k}*Se;
            Rotor{i}.C = Rotor{i}.C + Se'*Rotor{i}.Shaft{j}.R{k}'*Rotor{i}.Shaft{j}.Ke{k}*Rotor{i}.Shaft{j}.R{k}*Se * Rotor{i}.Shaft{j}.material.eta;
            
            %inertia
            Rotor{i}.M = Rotor{i}.M + Se'*Rotor{i}.Shaft{j}.R{k}'*Rotor{i}.Shaft{j}.Me{k}*Rotor{i}.Shaft{j}.R{k}*Se;
            Rotor{i}.G = Rotor{i}.G + Se'*Rotor{i}.Shaft{j}.R{k}'*Rotor{i}.Shaft{j}.Ge{k}*Rotor{i}.Shaft{j}.R{k}*Se * Rotor{i}.Speed;
            
            Rotor{i}.Fg = Rotor{i}.Fg + Se'*Rotor{i}.Shaft{j}.me(k)/2*repmat(info.gvec,2,1);
            
            Rotor{i}.F0 = Rotor{i}.F0 + Se'*Rotor{i}.Shaft{j}.R{k}'*Rotor{i}.Shaft{j}.F0{k};
        end
    end
    
    %now deal with the discs
    for j = 1:length(Rotor{i}.Disc)
        %Mesh.Ceate rotor->disc mapping matrix
        SDisc = IMapRotor(Rotor{i}.Disc{j}.iLocal,:);
        Rotor{i}.Disc{j}.S = SDisc;
        
        %disc compliance
        for k = 1:Rotor{i}.Disc{j}.Nt
            for l = 1:(Rotor{i}.Disc{j}.Nr-1)
                Se = Rotor{i}.Disc{j}.Se{k,l} * SDisc;
                Re = Rotor{i}.Disc{j}.Re{k,l};
                %stiffness/damping
                Rotor{i}.K = Rotor{i}.K + Se'*Re'*Rotor{i}.Disc{j}.Ke{k,l}*Re*Se;
                Rotor{i}.C = Rotor{i}.C + Se'*Re'*Rotor{i}.Disc{j}.Ke{k,l}*Re*Se * Rotor{i}.Disc{j}.material.eta;
                
                %only intertia terms
                Rotor{i}.M = Rotor{i}.M + Se'*Re'*Rotor{i}.Disc{j}.Me{k,l}*Re*Se;
                Rotor{i}.G = Rotor{i}.G + Se'*Re'*Rotor{i}.Disc{j}.Ge{k,l}*Re*Se * Rotor{i}.Speed;
                
                Rotor{i}.F0 = Rotor{i}.F0 + Se'*Re'*Rotor{i}.Disc{j}.Fe0{k,l};
            end
        end
        
        %edge stiffness
        SHub = Rotor{i}.Disc{j}.SHub*SDisc;
        for k = 1:Rotor{i}.Disc{j}.Nt
            SEdge = (Rotor{i}.Disc{j}.RHub{k,1}*SHub - Rotor{i}.Disc{j}.SNode{k,1}*SDisc);
            Rotor{i}.K = Rotor{i}.K + SEdge'*Rotor{i}.Disc{j}.KEdge*SEdge;
            Rotor{i}.C = Rotor{i}.C + SEdge'*Rotor{i}.Disc{j}.CEdge*SEdge;
        end
        
        %lump inertia at hub
        Rotor{i}.Fg = Rotor{i}.Fg + SHub'*Rotor{i}.Disc{j}.MHub*info.gvec;
        Rotor{i}.M  = Rotor{i}.M  + SHub'*Rotor{i}.Disc{j}.MHub*SHub;
        Rotor{i}.G  = Rotor{i}.G  + SHub'*Rotor{i}.Disc{j}.GHub*SHub;
        
        %root stiffness
        SRoot = Rotor{i}.Disc{j}.SRoot*SDisc;
        Rotor{i}.K  = Rotor{i}.K  + (SHub-SRoot)'*Rotor{i}.Disc{j}.KRoot*(SHub-SRoot);
        Rotor{i}.C  = Rotor{i}.C  + (SHub-SRoot)'*Rotor{i}.Disc{j}.CRoot*(SHub-SRoot);
        Rotor{i}.F0 = Rotor{i}.F0 + (SHub-SRoot)'*Rotor{i}.Disc{j}.F0Root;
    end
    
    Mesh.M  = Mesh.M  + Rotor{i}.S'*Rotor{i}.M*Rotor{i}.S;
    Mesh.G  = Mesh.G  + Rotor{i}.S'*Rotor{i}.G*Rotor{i}.S*Rotor{i}.Speed;
    Mesh.C  = Mesh.C  + Rotor{i}.S'*Rotor{i}.C*Rotor{i}.S;
    Mesh.K  = Mesh.K  + Rotor{i}.S'*Rotor{i}.K*Rotor{i}.S;
    Mesh.Fg  = Mesh.Fg + Rotor{i}.S'*Rotor{i}.Fg;
    Mesh.F0 = Mesh.F0 + Rotor{i}.S'*Rotor{i}.F0;
end

function [Bearing,Mesh] = setupbearings(Bearing,Mesh,Rotor,info)

NBearings = length(Bearing);

NRotorTot = 2*NBearings*info.NDofe;

Mesh.K = zeros(info.NDof);
Mesh.C = zeros(info.NDof);
Mesh.M = zeros(info.NDof);
Mesh.Fg = zeros(info.NDof,1);
Mesh.F0 = zeros(info.NDof,1);
Mesh.xInt = zeros(Mesh.NDofInt,1);
Mesh.u0 = zeros(2*NRotorTot,1);

IMapInput = eye(2*NRotorTot);
IMapInternal = eye(Mesh.NDofInt);

for i = 1:NBearings
    %mass of housing between bearing sets
    SBear = info.IMap(Bearing{i}.iGlobal,:);
    Mesh.M = Mesh.M + SBear'*Bearing{i}.M*SBear;
    Mesh.Fg = Mesh.Fg + SBear'*Bearing{i}.M*info.gvec;
    
    Bearing{i}.Sb = SBear; 
    
    %work out mapping matrices
    for j = 1:2 %out -> in
        if ~isnan(Bearing{i}.iRotor(j))
            iGlobal = Rotor{Bearing{i}.iRotor(j)}.iGlobal;
            iNode = Bearing{i}.iNode(j);
            Srot{j} = info.IMap(iGlobal((iNode-1)*info.NDofe+(1:info.NDofe)),:);
        else
            Srot{j} = zeros(info.NDofe,info.NDof);
        end
        
        Bearing{i}.Ui{j} = IMapInput(Bearing{i}.iInputi{j},:);
        Bearing{i}.Uo{j} = IMapInput(Bearing{i}.iInputo{j},:);
        Bearing{i}.V{j}  = IMapInternal(Bearing{i}.iInternal{j},:);
    end
    
    %outer side of bearing mass
    Bearing{i}.So{1} = Srot{1};
    Bearing{i}.Si{1} = SBear;
    
    %inner side of bearing mass
    Bearing{i}.So{2} = SBear;
    Bearing{i}.Si{2} = Srot{2};
    
    %now compute the stiffness matrices
    for j = 1:2
        Bearing{i}.S{j} = [Bearing{i}.Si{j}; Bearing{i}.So{j}];
        
        Mesh.K = Mesh.K + Bearing{i}.S{j}'*Bearing{i}.R{j}'*Bearing{i}.Kb{j}*Bearing{i}.R{j}*Bearing{i}.S{j};
        Mesh.C = Mesh.C + Bearing{i}.S{j}'*Bearing{i}.R{j}'*Bearing{i}.Cb{j}*Bearing{i}.R{j}*Bearing{i}.S{j};    
        Mesh.M = Mesh.M + Bearing{i}.S{j}'*Bearing{i}.R{j}'*Bearing{i}.Mb{j}*Bearing{i}.R{j}*Bearing{i}.S{j};    
        
        Mesh.F0 = Mesh.F0 + Bearing{i}.S{j}'*Bearing{i}.R{j}'*Bearing{i}.Fb{j};
        
        if ~isempty(Bearing{i}.xInt{j})
            Mesh.xInt = Mesh.xInt + Bearing{i}.V{j}'*Bearing{i}.xInt{j};
        end
        
        Mesh.u0 = Mesh.u0 + Bearing{i}.Ui{j}'*Bearing{i}.Ri{j}'*Bearing{i}.ui{j};
        Mesh.u0 = Mesh.u0 + Bearing{i}.Uo{j}'*Bearing{i}.Ro{j}'*Bearing{i}.uo{j};
    end
    
    %and finally work out the boundary nodes of the rotors
    for j = 1:2
        iRotor = Bearing{i}.iRotor(j);
        iNode  = Bearing{i}.iNode(j);
        if ~isnan(iRotor)
            Rotor{iRotor}.Bearing{end+1}.iBearing = i;
            Rotor{iRotor}.Bearing{end}.iNode = iNode;
            KFixed = Bearing{i}.R{j}'*Bearing{i}.Kb{j}*Bearing{i}.R{j};
            Rotor{iRotor}.Bearing{end}.iFixed = [];
            for k = 1:4
                if any(abs(KFixed(k,:))>0)
                    Rotor{iRotor}.Bearing{end}.iFixed(end+1,1) = k;
                end
            end
        end
    end
end

function [Excite,Mesh] = setupexcitation(Excite,Mesh,Rotor,Bearing,info)

%% Excitations
IMapExcite = eye(Mesh.NInput);

Mub = zeros(info.NDof);
Cub = zeros(info.NDof);
Kub = zeros(info.NDof);
Sub = zeros(info.NDof,Mesh.NInput);
uub = zeros(Mesh.NInput,1);

NBearingTot = length(Bearing)*info.NDofe;

Kgd = zeros(info.NDof,2*info.NDofBearing);
Cgd = zeros(info.NDof,2*info.NDofBearing);
Mgd = zeros(info.NDof,2*info.NDofBearing);
Sgd = zeros(2*NBearingTot,Mesh.NInput);
ugd = zeros(Mesh.NInput,1);

for i = 1:length(Excite)
    switch Excite{i}.Name
        case 'unbalance'
            iRotor = Excite{i}.iRotor;
            iDisc  = Excite{i}.iDisc;

            Sd = Rotor{iRotor}.Disc{iDisc}.SHub*Rotor{iRotor}.Disc{iDisc}.S*Rotor{iRotor}.S;
           
            Excite{i}.S = IMapExcite(Excite{i}.iExcite,:);
            
            Kub = Kub + Sd'*Excite{i}.K*Sd;
            Cub = Cub + Sd'*Excite{i}.C*Sd;
            Mub = Mub + Sd'*Excite{i}.M*Sd;
            
            Sub = Sub + Sd'*Excite{i}.S;
            
            uub = uub + Excite{i}.S'*Excite{i}.u;
        case 'Mesh.Gound'
            iBearing = Excite{i}.iBearing;
            
            %TODO: this needs updating to add deflection to inner/outer race
            Ub = Bearing{iBearing}.U;
            
            Sgnd = IMapExcite(Excite{i}.iExcite,:);
            
            for j = 1:2
                Excite{i}.U{j} = Sgnd((j-1)*info.NDofe + (1:info.NDofe),:);
                
                Sgd = Sgd + Ub{j}'*Excite{i}.U{j};

                Kgd = Kgd + Bearing{iBearing}.S{j}'*Bearing{i}.R{j}'*Bearing{iBearing}.Kb{j}*Bearing{i}.R{j}*Ub{j};
                Cgd = Cgd + Bearing{iBearing}.S{j}'*Bearing{i}.R{j}'*Bearing{iBearing}.Cb{j}*Bearing{i}.R{j}*Ub{j};
                Mgd = Mgd + Bearing{iBearing}.S{j}'*Bearing{i}.R{j}'*Bearing{iBearing}.Mb{j}*Bearing{i}.R{j}*Ub{j};
                
                ugd = ugd + Excite{i}.U'*Excite{i}.u;
            end
    end
end

Mesh.Kub = Kub;
Mesh.Cub = Cub;
Mesh.Mub = Mub;

Mesh.Kgd = Kgd;
Mesh.Cgd = Cgd;
Mesh.Mgd = Mgd;

Mesh.ugd = ugd;
Mesh.uub = uub;

Mesh.Sgd = Sgd;
Mesh.Sub = Sub;

Mesh.Ke = Mesh.Kub*Mesh.Sub + Mesh.Kgd*Mesh.Sgd;
Mesh.Ce = Mesh.Cub*Mesh.Sub + Mesh.Cgd*Mesh.Sgd;
Mesh.Me = Mesh.Mub*Mesh.Sub + Mesh.Mgd*Mesh.Sgd;
Mesh.ue = uub + ugd;
