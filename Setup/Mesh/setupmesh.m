function P = setupmesh(P)
NDofe = 4; %4 dof per elem

%find the number of rotor dof
NRotor = length(P.Rotor);
iDofCount = 0;
for i = 1:NRotor
    P.Rotor{i}.Bearing = {};
    NDofRotor(i) = P.Rotor{i}.NDof;
    P.Rotor{i}.iGlobal = iDofCount + (1:P.Rotor{i}.NDof);
    iDofCount = iDofCount + P.Rotor{i}.NDof;
end
NDofRotorTot = sum(NDofRotor);

%and the bearings
NBearings = length(P.Bearing);
for i = 1:NBearings
   P.Bearing{i}.iGlobal = iDofCount + (1:NDofe);
   iDofCount = iDofCount + NDofe;
end
NDofBearing = NDofe * NBearings;

NDofTot = NDofRotorTot + NDofBearing;
IMapGlobal = eye(NDofTot);

%% Begin with the rotors
Kr = zeros(NDofTot);
Cr = zeros(NDofTot);
Mr = zeros(NDofTot);
Gr = zeros(NDofTot);
F0r = zeros(NDofTot,1);
Fgr = zeros(NDofTot,1);

for i = 1:NRotor
    %create the global -> rotor mapping matrix
    Sr  = IMapGlobal(P.Rotor{i}.iGlobal,:);
    P.Rotor{i}.S = Sr;
        
    P.Rotor{i}.M = zeros(P.Rotor{i}.NDof);
    P.Rotor{i}.G = zeros(P.Rotor{i}.NDof);
    P.Rotor{i}.C = zeros(P.Rotor{i}.NDof);
    P.Rotor{i}.K = zeros(P.Rotor{i}.NDof);
    P.Rotor{i}.Fg = zeros(P.Rotor{i}.NDof,1);
    P.Rotor{i}.F0 = zeros(P.Rotor{i}.NDof,1);
    
    IMapRotor = eye(P.Rotor{i}.NDof);
    
    %deal with the shafts first
    for j = 1:length(P.Rotor{i}.Shaft)
        SShaft = IMapRotor(P.Rotor{i}.Shaft{j}.iLocal,:);
        P.Rotor{i}.Shaft{j}.S = SShaft;
        
        for k = 1:(P.Rotor{i}.Shaft{j}.Nz-1)
            %rotor->shaft element mapping matrix
            Se = P.Rotor{i}.Shaft{j}.Se{k}*SShaft;
                        
            %stiffness/damping
            P.Rotor{i}.K = P.Rotor{i}.K + Se'*P.Rotor{i}.Shaft{j}.R{k}'*P.Rotor{i}.Shaft{j}.Ke{k}*P.Rotor{i}.Shaft{j}.R{k}*Se;
            P.Rotor{i}.C = P.Rotor{i}.C + Se'*P.Rotor{i}.Shaft{j}.R{k}'*P.Rotor{i}.Shaft{j}.Ke{k}*P.Rotor{i}.Shaft{j}.R{k}*Se * P.Rotor{i}.Shaft{j}.Material.eta;
            
            %inertia
            P.Rotor{i}.M = P.Rotor{i}.M + Se'*P.Rotor{i}.Shaft{j}.R{k}'*P.Rotor{i}.Shaft{j}.Me{k}*P.Rotor{i}.Shaft{j}.R{k}*Se;
            P.Rotor{i}.G = P.Rotor{i}.G + Se'*P.Rotor{i}.Shaft{j}.R{k}'*P.Rotor{i}.Shaft{j}.Ge{k}*P.Rotor{i}.Shaft{j}.R{k}*Se * P.Rotor{i}.Speed;
            
            P.Rotor{i}.Fg = P.Rotor{i}.Fg + Se'*P.Rotor{i}.Shaft{j}.me(k)/2*repmat([P.g; 0; 0],2,1);
            
            P.Rotor{i}.F0 = P.Rotor{i}.F0 + Se'*P.Rotor{i}.Shaft{j}.R{k}'*P.Rotor{i}.Shaft{j}.F0{k};
        end
    end
    
    %now deal with the discs
    for j = 1:length(P.Rotor{i}.Disc)
        %create rotor->disc mapping matrix
        SDisc = IMapRotor(P.Rotor{i}.Disc{j}.iLocal,:);
        P.Rotor{i}.Disc{j}.S = SDisc;
        
        %disc compliance
        for k = 1:P.Rotor{i}.Disc{j}.Nt
            for l = 1:(P.Rotor{i}.Disc{j}.Nr-1)
                Se = P.Rotor{i}.Disc{j}.Se{k,l} * SDisc;
                Re = P.Rotor{i}.Disc{j}.Re{k,l};
                %stiffness/damping
                P.Rotor{i}.K = P.Rotor{i}.K + Se'*Re'*P.Rotor{i}.Disc{j}.Ke{k,l}*Re*Se;
                P.Rotor{i}.C = P.Rotor{i}.C + Se'*Re'*P.Rotor{i}.Disc{j}.Ke{k,l}*Re*Se * P.Rotor{i}.Disc{j}.Material.eta;
                
                %only intertia terms
                P.Rotor{i}.M = P.Rotor{i}.M + Se'*Re'*P.Rotor{i}.Disc{j}.Me{k,l}*Re*Se;
                P.Rotor{i}.G = P.Rotor{i}.G + Se'*Re'*P.Rotor{i}.Disc{j}.Ge{k,l}*Re*Se * P.Rotor{i}.Speed;
                
                P.Rotor{i}.F0 = P.Rotor{i}.F0 + Se'*Re'*P.Rotor{i}.Disc{j}.Fe0{k,l};
            end
        end
        
        %edge stiffness
        SHub = P.Rotor{i}.Disc{j}.SHub*SDisc;
        for k = 1:P.Rotor{i}.Disc{j}.Nt
            SEdge = (P.Rotor{i}.Disc{j}.RHub{k,1}*SHub - P.Rotor{i}.Disc{j}.SNode{k,1}*SDisc);
            P.Rotor{i}.K = P.Rotor{i}.K + SEdge'*P.Rotor{i}.Disc{j}.KEdge*SEdge;
            P.Rotor{i}.C = P.Rotor{i}.C + SEdge'*P.Rotor{i}.Disc{j}.CEdge*SEdge;
        end
        
        %lump inertia at hub
        P.Rotor{i}.Fg = P.Rotor{i}.Fg + SHub'*P.Rotor{i}.Disc{j}.MHub*[P.g;0;0];
        P.Rotor{i}.M  = P.Rotor{i}.M  + SHub'*P.Rotor{i}.Disc{j}.MHub*SHub;
        P.Rotor{i}.G  = P.Rotor{i}.G  + SHub'*P.Rotor{i}.Disc{j}.GHub*SHub;
        
        %root stiffness
        SRoot = P.Rotor{i}.Disc{j}.SRoot*SDisc;
        P.Rotor{i}.K  = P.Rotor{i}.K  + (SHub-SRoot)'*P.Rotor{i}.Disc{j}.KRoot*(SHub-SRoot);
        P.Rotor{i}.C  = P.Rotor{i}.C  + (SHub-SRoot)'*P.Rotor{i}.Disc{j}.CRoot*(SHub-SRoot);
        P.Rotor{i}.F0 = P.Rotor{i}.F0 + (SHub-SRoot)'*P.Rotor{i}.Disc{j}.F0Root;
    end
    
    Mr  = Mr  + P.Rotor{i}.S'*P.Rotor{i}.M*P.Rotor{i}.S;
    Gr  = Gr  + P.Rotor{i}.S'*P.Rotor{i}.G*P.Rotor{i}.S*P.Rotor{i}.Speed;
    Cr  = Cr  + P.Rotor{i}.S'*P.Rotor{i}.C*P.Rotor{i}.S;
    Kr  = Kr  + P.Rotor{i}.S'*P.Rotor{i}.K*P.Rotor{i}.S;
    Fgr = Fgr + P.Rotor{i}.S'*P.Rotor{i}.Fg;
    F0r = F0r + P.Rotor{i}.S'*P.Rotor{i}.F0;
end

P.Mesh.Rotor.M = Mr;
P.Mesh.Rotor.G = Gr;
P.Mesh.Rotor.C = Cr;
P.Mesh.Rotor.K = Kr;
P.Mesh.Rotor.Fg = Fgr;
P.Mesh.Rotor.F0 = F0r;

%% Move onto the bearings

%work out the number of internal states
NInternal = zeros(NBearings,2);
for i = 1:NBearings
    NInternal(i,:) = P.Bearing{i}.NDofInt;
end
NInternalTot = sum(NInternal(:));
NBearingTot = 2*NBearings*NDofe;

Kb = zeros(NDofTot);
Cb = zeros(NDofTot);
Mb = zeros(NDofTot);
Fgb = zeros(NDofTot,1);
Fb = zeros(NDofTot,1);
xInt = zeros(NInternalTot,1);
u0 = zeros(2*NBearingTot,1);

IMapInput = eye(2*NBearingTot);
IMapInternal = eye(NInternalTot);

for i = 1:NBearings
    %mass of housing between bearing sets
    SBear = IMapGlobal(P.Bearing{i}.iGlobal,:);
    Mb = Mb + SBear'*P.Bearing{i}.M*SBear;
    Fgb = Fgb + SBear'*P.Bearing{i}.M*[P.g; 0; 0];
    
    P.Bearing{i}.Sb = SBear; 
    
    %work out mapping matrices
    for j = 1:2 %out -> in
        if ~isnan(P.Bearing{i}.iRotor(j))
            iGlobal = P.Rotor{P.Bearing{i}.iRotor(j)}.iGlobal;
            iNode = P.Bearing{i}.iNode(j);
            Srot{j} = IMapGlobal(iGlobal((iNode-1)*NDofe+(1:NDofe)),:);
        else
            Srot{j} = zeros(NDofe,NDofTot);
        end
        
        P.Bearing{i}.Ui{j} = IMapInput(P.Bearing{i}.iInputi{j},:);
        P.Bearing{i}.Uo{j} = IMapInput(P.Bearing{i}.iInputo{j},:);
        P.Bearing{i}.V{j}  = IMapInternal(P.Bearing{i}.iInternal{j},:);
    end
    
    %outer side of bearing mass
    P.Bearing{i}.So{1} = Srot{1};
    P.Bearing{i}.Si{1} = SBear;
    
    %inner side of bearing mass
    P.Bearing{i}.So{2} = SBear;
    P.Bearing{i}.Si{2} = Srot{2};
    
    %now compute the stiffness matrices
    for j = 1:2
        P.Bearing{i}.S{j} = [P.Bearing{i}.Si{j}; P.Bearing{i}.So{j}];
        
        Kb = Kb + P.Bearing{i}.S{j}'*P.Bearing{i}.R{j}'*P.Bearing{i}.Kb{j}*P.Bearing{i}.R{j}*P.Bearing{i}.S{j};
        Cb = Cb + P.Bearing{i}.S{j}'*P.Bearing{i}.R{j}'*P.Bearing{i}.Cb{j}*P.Bearing{i}.R{j}*P.Bearing{i}.S{j};    
        Mb = Mb + P.Bearing{i}.S{j}'*P.Bearing{i}.R{j}'*P.Bearing{i}.Mb{j}*P.Bearing{i}.R{j}*P.Bearing{i}.S{j};    
        
        Fb = Fb + P.Bearing{i}.S{j}'*P.Bearing{i}.R{j}'*P.Bearing{i}.Fb{j};
        
        if ~isempty(P.Bearing{i}.xInt{j})
            xInt = xInt + P.Bearing{i}.V{j}'*P.Bearing{i}.xInt{j};
        end
        
        u0 = u0 + P.Bearing{i}.Ui{j}'*P.Bearing{i}.Ri{j}'*P.Bearing{i}.ui{j};
        u0 = u0 + P.Bearing{i}.Uo{j}'*P.Bearing{i}.Ro{j}'*P.Bearing{i}.uo{j};
    end
    
    %and finally work out the boundary nodes of the rotors
    for j = 1:2
        iRotor = P.Bearing{i}.iRotor(j);
        iNode  = P.Bearing{i}.iNode(j);
        if ~isnan(iRotor)
            P.Rotor{iRotor}.Bearing{end+1}.iBearing = i;
            P.Rotor{iRotor}.Bearing{end}.iNode = iNode;
            P.Rotor{iRotor}.Bearing{end}.iActive = findrows(P.Bearing{i}.Ri{j}(P.Bearing{i}.bActive{j},:));
        end
    end
end

P.Mesh.Bearing.K = Kb;
P.Mesh.Bearing.C = Cb;
P.Mesh.Bearing.M = Mb;
P.Mesh.Bearing.Fg = Fgb;
P.Mesh.Bearing.F0 = Fb;
P.Mesh.Bearing.u0 = u0;
P.Mesh.Bearing.xInt = xInt;

%% Excitations
for i = 1:length(P.Excite)
    NExcite(i) = P.Excite{i}.NInput;
end
NExciteTot = sum(NExcite(:));

IMapExcite = eye(NExciteTot);

Mub = zeros(NDofTot);
Cub = zeros(NDofTot);
Kub = zeros(NDofTot);
Sub = zeros(NDofTot,NExciteTot);
uub = zeros(NExciteTot,1);

Kgd = zeros(NDofTot,2*NBearingTot);
Cgd = zeros(NDofTot,2*NBearingTot);
Mgd = zeros(NDofTot,2*NBearingTot);
Sgd = zeros(2*NBearingTot,NExciteTot);
ugd = zeros(NExciteTot,1);

for i = 1:length(P.Excite)
    switch P.Excite{i}.Name
        case 'unbalance'
            iRotor = P.Excite{i}.iRotor;
            iDisc  = P.Excite{i}.iDisc;

            Sd = P.Rotor{iRotor}.Disc{iDisc}.SHub*P.Rotor{iRotor}.Disc{iDisc}.S*P.Rotor{iRotor}.S;
           
            P.Excite{i}.S = IMapExcite(P.Excite{i}.iExcite,:);
            
            Kub = Kub + Sd'*P.Excite{i}.K*Sd;
            Cub = Cub + Sd'*P.Excite{i}.C*Sd;
            Mub = Mub + Sd'*P.Excite{i}.M*Sd;
            
            Sub = Sub + Sd'*P.Excite{i}.S;
            
            uub = uub + P.Excite{i}.S'*P.Excite{i}.u;
        case 'ground'
            iBearing = P.Excite{i}.iBearing;
            
            %TODO: this needs updating to add deflection to inner/outer race
            Ub = P.Bearing{iBearing}.U;
            
            Sgnd = IMapExcite(P.Excite{i}.iExcite,:);
            
            for j = 1:2
                P.Excite{i}.U{j} = Sgnd((j-1)*NDofe + (1:NDofe),:);
                
                Sgd = Sgd + Ub{j}'*P.Excite{i}.U{j};

                Kgd = Kgd + P.Bearing{iBearing}.S{j}'*P.Bearing{i}.R{j}'*P.Bearing{iBearing}.Kb{j}*P.Bearing{i}.R{j}*Ub{j};
                Cgd = Cgd + P.Bearing{iBearing}.S{j}'*P.Bearing{i}.R{j}'*P.Bearing{iBearing}.Cb{j}*P.Bearing{i}.R{j}*Ub{j};
                Mgd = Mgd + P.Bearing{iBearing}.S{j}'*P.Bearing{i}.R{j}'*P.Bearing{iBearing}.Mb{j}*P.Bearing{i}.R{j}*Ub{j};
                
                ugd = ugd + P.Excite{i}.U'*P.Excite{i}.u;
            end
    end
end

P.Mesh.Excite.Kub = Kub;
P.Mesh.Excite.Cub = Cub;
P.Mesh.Excite.Mub = Mub;

P.Mesh.Excite.Kgd = Kgd;
P.Mesh.Excite.Cgd = Cgd;
P.Mesh.Excite.Mgd = Mgd;

P.Mesh.Excite.ugd = ugd;
P.Mesh.Excite.uub = uub;

P.Mesh.Excite.Sgd = Sgd;
P.Mesh.Excite.Sub = Sub;

P.Mesh.Excite.Ke = P.Mesh.Excite.Kub*P.Mesh.Excite.Sub + P.Mesh.Excite.Kgd*P.Mesh.Excite.Sgd;
P.Mesh.Excite.Ce = P.Mesh.Excite.Cub*P.Mesh.Excite.Sub + P.Mesh.Excite.Cgd*P.Mesh.Excite.Sgd;
P.Mesh.Excite.Me = P.Mesh.Excite.Mub*P.Mesh.Excite.Sub + P.Mesh.Excite.Mgd*P.Mesh.Excite.Sgd;
P.Mesh.Excite.ue = uub + ugd;

%% Combined
P.Mesh.M  = P.Mesh.Rotor.M  + P.Mesh.Bearing.M;
P.Mesh.G  = P.Mesh.Rotor.G;                    
P.Mesh.K  = P.Mesh.Rotor.K  + P.Mesh.Bearing.K;
P.Mesh.C  = P.Mesh.Rotor.C  + P.Mesh.Bearing.C;
P.Mesh.Fg = P.Mesh.Rotor.Fg + P.Mesh.Bearing.Fg;
P.Mesh.F0 = P.Mesh.Rotor.F0 + P.Mesh.Bearing.F0;
P.Mesh.A  = eye(NDofTot);

%% Store some useful numbers
P.Mesh.NDofInt = NInternalTot;
P.Mesh.NDof = NDofTot;
P.Mesh.NDofTot = NDofTot+NInternalTot;
P.Mesh.NExcite = NExciteTot;

P.Mesh.Rotor.NDof = NDofRotor;

P.Mesh.Bearing.NDof    = NDofBearing;
P.Mesh.Bearing.NDofInt = NInternalTot;