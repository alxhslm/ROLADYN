function P = setupmesh(P)
NDofe = 4; %4 dof per elem

%find the number of rotor dof
NRotor = length(P.Rotor);
iDofCount = 0;
NDofRotor = zeros(NRotor,1);
for i = 1:NRotor
    P.Rotor{i}.Bearing = {};
    NDofRotor(i) = P.Rotor{i}.NDof;
    P.Rotor{i}.iGlobal = iDofCount + (1:P.Rotor{i}.NDof);
    iDofCount = iDofCount + P.Rotor{i}.NDof;
end
NDofRotorTot = sum(NDofRotor);

%and the stator
NStator = length(P.Stator);
NDofStator = zeros(NStator,1);
for i = 1:NStator
   P.Stator{i}.Bearing = {};
   NDofStator(i) = P.Stator{i}.NDof;
   P.Stator{i}.iGlobal = iDofCount + (1:P.Stator{i}.NDof);
   iDofCount = iDofCount + P.Stator{i}.NDof;
end
NDofStatorTot = sum(NDofStator);

NDofTot = NDofRotorTot + NDofStatorTot;
IMapGlobal = eye(NDofTot);

%% Begin with the rotors
Kr = zeros(NDofTot);
Cr = zeros(NDofTot);
Mr = zeros(NDofTot);
Gr = zeros(NDofTot);
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
    
    IMapRotor = eye(P.Rotor{i}.NDof);
    
    %deal with the shafts first
    for j = 1:length(P.Rotor{i}.Shaft)
        SShaft = IMapRotor(P.Rotor{i}.Shaft{j}.iLocal,:);
        P.Rotor{i}.Shaft{j}.S = SShaft;
        
        for k = 1:(P.Rotor{i}.Shaft{j}.Mesh.Nz-1)
            %rotor->shaft element mapping matrix
            Se = P.Rotor{i}.Shaft{j}.Element{k}.S*SShaft;
            Re = P.Rotor{i}.Shaft{j}.Element{k}.R;
                        
            %stiffness/damping
            P.Rotor{i}.K = P.Rotor{i}.K + Se'*Re'*P.Rotor{i}.Shaft{j}.Element{k}.K*Re*Se;
            P.Rotor{i}.C = P.Rotor{i}.C + Se'*Re'*P.Rotor{i}.Shaft{j}.Element{k}.C*Re*Se;
            
            %inertia
            P.Rotor{i}.M = P.Rotor{i}.M + Se'*Re'*P.Rotor{i}.Shaft{j}.Element{k}.M*Re*Se;
            P.Rotor{i}.G = P.Rotor{i}.G + Se'*Re'*P.Rotor{i}.Shaft{j}.Element{k}.G*Re*Se * P.Rotor{i}.Speed;
            
            P.Rotor{i}.Fg = P.Rotor{i}.Fg + Se'*Re'*P.Rotor{i}.Shaft{j}.Element{k}.M*Re*repmat([P.g; 0; 0],2,1);
        end
    end
    
    %now deal with the discs
    for j = 1:length(P.Rotor{i}.Disc)
        %create rotor->disc mapping matrix
        SDisc = IMapRotor(P.Rotor{i}.Disc{j}.iLocal,:);
        P.Rotor{i}.Disc{j}.S = SDisc;
        
        SHub = P.Rotor{i}.Disc{j}.Hub.S*SDisc;
        SRoot = P.Rotor{i}.Disc{j}.Root.S*SDisc;
        
        if strcmp(P.Rotor{i}.Disc{j}.Type,'Flexible')
            %disc compliance
            for k = 1:P.Rotor{i}.Disc{j}.Mesh.Nt
                for l = 1:(P.Rotor{i}.Disc{j}.Mesh.Nr-1)
                    Se = P.Rotor{i}.Disc{j}.Element{k,l}.S * SDisc;
                    Re = P.Rotor{i}.Disc{j}.Element{k,l}.R;
                    %stiffness/damping
                    P.Rotor{i}.K = P.Rotor{i}.K + Se'*Re'*P.Rotor{i}.Disc{j}.Element{k,l}.K*Re*Se;
                    P.Rotor{i}.C = P.Rotor{i}.C + Se'*Re'*P.Rotor{i}.Disc{j}.Element{k,l}.C*Re*Se;
                    
                    %only intertia terms
                    P.Rotor{i}.M = P.Rotor{i}.M + Se'*Re'*P.Rotor{i}.Disc{j}.Element{k,l}.M*Re*Se;
                    P.Rotor{i}.G = P.Rotor{i}.G + Se'*Re'*P.Rotor{i}.Disc{j}.Element{k,l}.G*Re*Se * P.Rotor{i}.Speed;
                end
            end
            
            %edge stiffness
            for k = 1:P.Rotor{i}.Disc{j}.Mesh.Nt
                SEdge = (P.Rotor{i}.Disc{j}.Mesh.RHub{k,1}*SHub - P.Rotor{i}.Disc{j}.Mesh.SNode{k,1}*SDisc);
                P.Rotor{i}.K = P.Rotor{i}.K + SEdge'*P.Rotor{i}.Disc{j}.Edge.K*SEdge;
                P.Rotor{i}.C = P.Rotor{i}.C + SEdge'*P.Rotor{i}.Disc{j}.Edge.C*SEdge;
            end
        end
        
        %lump inertia at hub
        P.Rotor{i}.Fg = P.Rotor{i}.Fg + SHub'*(P.Rotor{i}.Disc{j}.Hub.M + P.Rotor{i}.Disc{j}.Ring.M)*[P.g;0;0];
        P.Rotor{i}.M  = P.Rotor{i}.M  + SHub'*(P.Rotor{i}.Disc{j}.Hub.M + P.Rotor{i}.Disc{j}.Ring.M)*SHub;
        P.Rotor{i}.G  = P.Rotor{i}.G  + SHub'*(P.Rotor{i}.Disc{j}.Hub.G + P.Rotor{i}.Disc{j}.Ring.G)*SHub;
        
        %root stiffness
        P.Rotor{i}.K  = P.Rotor{i}.K  + (SHub-SRoot)'*P.Rotor{i}.Disc{j}.Root.K*(SHub-SRoot);
        P.Rotor{i}.C  = P.Rotor{i}.C  + (SHub-SRoot)'*P.Rotor{i}.Disc{j}.Root.C*(SHub-SRoot);
    end
    
    Mr  = Mr  + P.Rotor{i}.S'*P.Rotor{i}.M*P.Rotor{i}.S;
    Gr  = Gr  + P.Rotor{i}.S'*P.Rotor{i}.G*P.Rotor{i}.S*P.Rotor{i}.Speed;
    Cr  = Cr  + P.Rotor{i}.S'*P.Rotor{i}.C*P.Rotor{i}.S;
    Kr  = Kr  + P.Rotor{i}.S'*P.Rotor{i}.K*P.Rotor{i}.S;
    Fgr = Fgr + P.Rotor{i}.S'*P.Rotor{i}.Fg;
end

P.Mesh.Rotor.M = Mr;
P.Mesh.Rotor.G = Gr;
P.Mesh.Rotor.C = Cr;
P.Mesh.Rotor.K = Kr;
P.Mesh.Rotor.Fg = Fgr;

%% And the stator
Ks = zeros(NDofTot);
Cs = zeros(NDofTot);
Ms = zeros(NDofTot);
Fgs = zeros(NDofTot,1);

for i = 1:NStator
    %create the global -> rotor mapping matrix
    Ss  = IMapGlobal(P.Stator{i}.iGlobal,:);
    P.Stator{i}.S = Ss;
    
    P.Stator{i}.Fg = P.Stator{i}.M*[P.g;zeros(P.Stator{i}.NDof-2,1)];
       
    Ms  = Ms  + P.Stator{i}.S'*P.Stator{i}.M*P.Stator{i}.S;
    Cs  = Cs  + P.Stator{i}.S'*P.Stator{i}.C*P.Stator{i}.S;
    Ks  = Ks  + P.Stator{i}.S'*P.Stator{i}.K*P.Stator{i}.S;
    Fgs = Fgs + P.Stator{i}.S'*P.Stator{i}.Fg;
end

P.Mesh.Stator.M = Ms;
P.Mesh.Stator.C = Cs;
P.Mesh.Stator.K = Ks;
P.Mesh.Stator.Fg = Fgs;

%% Move onto the bearings
NBearings = length(P.Bearing);

NInputBearing = 2*NDofe*NBearings;
Fb = zeros(NInputBearing,1);
Kb = zeros(NInputBearing);
Cb = zeros(NInputBearing);
Mb = zeros(NInputBearing);

Fb_lin = zeros(NInputBearing,1);
Kb_lin = zeros(NInputBearing);
Cb_lin = zeros(NInputBearing);
Mb_lin = zeros(NInputBearing);

IMapInput = eye(NInputBearing);

%NB. order of nodes: inner, outer 
Sb = zeros(0,NDofTot);
for i = 1:NBearings
    bGround = false;
    %work out mapping matrices
    for j = 1:2
        switch P.Bearing{i}.Node{j}.Type
            case 'ground'
                Sio{j} = zeros(NDofe,NDofTot);
                z(j) = NaN;
                bGround = true;
            case 'rotor'
                iNode = P.Bearing{i}.Node{j}.iNode;
                Sio{j} = P.Rotor{P.Bearing{i}.Node{j}.iRotor}.S((iNode-1)*NDofe+(1:NDofe),:);
                z(j) =  P.Rotor{P.Bearing{i}.Node{j}.iRotor}.Nodes(iNode);
            case 'stator'
                Sio{j} = P.Stator{P.Bearing{i}.Node{j}.iStator}.S(1:NDofe,:);
                z(j) =  P.Stator{P.Bearing{i}.Node{j}.iStator}.z;
        end
    end

    bLinear = strcmp(P.Bearing{i}.Model,'linear');
    
    Rio = {P.Bearing{i}.Ri,P.Bearing{i}.Ro};
    
    P.Bearing{i}.zi = z(1);
    P.Bearing{i}.zo = z(2);

    P.Bearing{i}.Ti = axial_offset(P.Bearing{i}.z-z(1));
    P.Bearing{i}.To = axial_offset(P.Bearing{i}.z-z(2));
        
    P.Bearing{i}.Si = P.Bearing{i}.Ti*Sio{1};
    P.Bearing{i}.So = P.Bearing{i}.To*Sio{2};
    P.Bearing{i}.S = [P.Bearing{i}.Si; P.Bearing{i}.So];
    Sb = [Sb; P.Bearing{i}.S];

    P.Bearing{i}.Ui = IMapInput(P.Bearing{i}.iInputi,:);
    P.Bearing{i}.Uo = IMapInput(P.Bearing{i}.iInputo,:);
    P.Bearing{i}.U = [P.Bearing{i}.Ui;P.Bearing{i}.Uo];
    P.Bearing{i}.R = blkdiag(P.Bearing{i}.Ri,P.Bearing{i}.Ro);

    %now compute the stiffness matrices  
    Fb = Fb + P.Bearing{i}.U'*P.Bearing{i}.R'*P.Bearing{i}.F0;
    Kb = Kb + P.Bearing{i}.U'*P.Bearing{i}.R'*P.Bearing{i}.Kb*P.Bearing{i}.R*P.Bearing{i}.U;
    Cb = Cb + P.Bearing{i}.U'*P.Bearing{i}.R'*P.Bearing{i}.Cb*P.Bearing{i}.R*P.Bearing{i}.U;
    Mb = Mb + P.Bearing{i}.U'*P.Bearing{i}.R'*P.Bearing{i}.Mb*P.Bearing{i}.R*P.Bearing{i}.U;
    
    if bLinear
        Fb_lin = Fb_lin + P.Bearing{i}.U'*P.Bearing{i}.R'*P.Bearing{i}.F0;
        Kb_lin = Kb_lin + P.Bearing{i}.U'*P.Bearing{i}.R'*P.Bearing{i}.Kb*P.Bearing{i}.R*P.Bearing{i}.U;
        Cb_lin = Cb_lin + P.Bearing{i}.U'*P.Bearing{i}.R'*P.Bearing{i}.Cb*P.Bearing{i}.R*P.Bearing{i}.U;
        Mb_lin = Mb_lin + P.Bearing{i}.U'*P.Bearing{i}.R'*P.Bearing{i}.Mb*P.Bearing{i}.R*P.Bearing{i}.U;
    end
    
    %and finally work out the boundary nodes of the rotors
    for j = 1:2
        switch P.Bearing{i}.Node{j}.Type
            case 'rotor'
                iRotor = P.Bearing{i}.Node{j}.iRotor;
                iNode  = P.Bearing{i}.Node{j}.iNode;
                
                P.Rotor{iRotor}.Bearing{end+1}.iBearing = i;
                P.Rotor{iRotor}.Bearing{end}.iNode = iNode;
                P.Rotor{iRotor}.Bearing{end}.iActive = findrows(Rio{j}(P.Bearing{i}.bActive,:));
                P.Rotor{iRotor}.Bearing{end}.bRigid = P.Bearing{i}.bRigid(P.Bearing{i}.bActive);
                P.Rotor{iRotor}.Bearing{end}.iNodeBearing = j;
                P.Rotor{iRotor}.Bearing{end}.bLinear = bLinear;
                P.Rotor{iRotor}.Bearing{end}.bGround = bGround;
            case 'stator'
                iStator = P.Bearing{i}.Node{j}.iStator;
                
                P.Stator{iStator}.Bearing{end+1}.iBearing = i;
                P.Stator{iStator}.Bearing{end}.iActive = findrows(Rio{j}(P.Bearing{i}.bActive,:));
                P.Stator{iStator}.Bearing{end}.bRigid = P.Bearing{i}.bRigid(P.Bearing{i}.bActive);
        end
    end
        
    P.Bearing{i}.bLinear = bLinear;
    P.Bearing{i}.bGround = bGround;
end

P.Mesh.Bearing.Fb = Fb;
P.Mesh.Bearing.Kb = Kb;
P.Mesh.Bearing.Cb = Cb;
P.Mesh.Bearing.Mb = Mb;
P.Mesh.Bearing.F  = Sb'*Fb;
P.Mesh.Bearing.K  = Sb'*Kb*Sb;
P.Mesh.Bearing.C  = Sb'*Cb*Sb;
P.Mesh.Bearing.M  = Sb'*Mb*Sb;

P.Mesh.Bearing.Lin.Fb = Fb_lin;
P.Mesh.Bearing.Lin.Kb = Kb_lin;
P.Mesh.Bearing.Lin.Cb = Cb_lin;
P.Mesh.Bearing.Lin.Mb = Mb_lin;
P.Mesh.Bearing.Lin.F  = Sb'*Fb_lin;
P.Mesh.Bearing.Lin.K  = Sb'*Kb_lin*Sb;
P.Mesh.Bearing.Lin.C  = Sb'*Cb_lin*Sb;
P.Mesh.Bearing.Lin.M  = Sb'*Mb_lin*Sb;

P.Mesh.Bearing.NL.Fb = Fb-Fb_lin;
P.Mesh.Bearing.NL.Kb = Kb-Kb_lin;
P.Mesh.Bearing.NL.Cb = Cb-Cb_lin;
P.Mesh.Bearing.NL.Mb = Mb-Mb_lin;
P.Mesh.Bearing.NL.F = Sb'*(Fb-Fb_lin);
P.Mesh.Bearing.NL.K = Sb'*(Kb-Kb_lin)*Sb;
P.Mesh.Bearing.NL.C = Sb'*(Cb-Cb_lin)*Sb;
P.Mesh.Bearing.NL.M = Sb'*(Mb-Mb_lin)*Sb;

P.Mesh.Bearing.S = Sb;

%% Excitations
NExcInput = zeros(length(P.Excite),1);
for i = 1:length(P.Excite)
    NExcInput(i) = P.Excite{i}.NInput;
end
NExcInputTot = sum(NExcInput(:));

IMapExcite = eye(NExcInputTot);

Me = zeros(NExcInputTot);
Ce = zeros(NExcInputTot);
Ke = zeros(NExcInputTot);
usync = zeros(NExcInputTot,1);
uasync = zeros(NExcInputTot,1);
Se = zeros(0,NDofTot);

for i = 1:NBearings
    P.Bearing{i}.NExcite = 0;
    P.Bearing{i}.Ue = zeros(0,NExcInputTot);
end

for i = 1:length(P.Excite)
    Ue = IMapExcite(P.Excite{i}.iExcite,:);
    P.Excite{i}.U = Ue;

    switch P.Excite{i}.Type
        case 'unbalance'
            iRotor = P.Excite{i}.iRotor;
            iDisc  = P.Excite{i}.iDisc;
            P.Excite{i}.S = P.Rotor{iRotor}.Disc{iDisc}.Hub.S(1:2,:)*P.Rotor{iRotor}.Disc{iDisc}.S*P.Rotor{iRotor}.S;
            Me = Me + Ue'*P.Excite{i}.M*Ue;
        case 'skew'
            iRotor = P.Excite{i}.iRotor;
            iDisc  = P.Excite{i}.iDisc;
            P.Excite{i}.S = P.Rotor{iRotor}.Disc{iDisc}.Hub.S(3:4,:)*P.Rotor{iRotor}.Disc{iDisc}.S*P.Rotor{iRotor}.S;
            Me = Me + Ue'*P.Excite{i}.M*Ue;
        case 'shaker'
            iStator = P.Excite{i}.iStator;
            P.Excite{i}.S = P.Stator{iStator}.S(1:2,:);
            Ke = Ke + Ue'*P.Excite{i}.K*Ue;
            Ce = Ce + Ue'*P.Excite{i}.C*Ue;
            Me = Me + Ue'*P.Excite{i}.M*Ue;
        case 'bearing'
            iBearing = P.Excite{i}.iBearing;
            P.Excite{i}.S = zeros(P.Excite{i}.NInput,NDofTot);           
            P.Bearing{iBearing}.NExcite = P.Excite{i}.NInput;
            P.Bearing{iBearing}.Ue = Ue;
    end
   
    Se = [Se; P.Excite{i}.S];
    switch P.Excite{i}.Mode
        case 'sync'
           usync = usync + Ue'*P.Excite{i}.u;
        case 'async'
           uasync = uasync + Ue'*P.Excite{i}.u;
        otherwise
            error('Unknown excitation mode "%s"',P.Excite{i}.Mode)
    end
end

P.Mesh.Excite.Ke = Ke;
P.Mesh.Excite.Ce = Ce;
P.Mesh.Excite.Me = Me;

P.Mesh.Excite.K = Se'*Ke;
P.Mesh.Excite.C = Se'*Ce;
P.Mesh.Excite.M = Se'*Me;
P.Mesh.Excite.S = Se;

P.Mesh.Excite.uSync  = usync;
P.Mesh.Excite.uAsync = uasync;

%% Combined
P.Mesh.M  = P.Mesh.Rotor.M + P.Mesh.Stator.M;
P.Mesh.G  = P.Mesh.Rotor.G;                    
P.Mesh.K  = P.Mesh.Rotor.K  + P.Mesh.Bearing.K + P.Mesh.Stator.K;
P.Mesh.C  = P.Mesh.Rotor.C  + P.Mesh.Bearing.C + P.Mesh.Stator.C;
P.Mesh.Fg = P.Mesh.Rotor.Fg + P.Mesh.Stator.Fg;
P.Mesh.A  = eye(NDofTot);

%% Store some useful numbers
P.Mesh.NDof = NDofTot;
P.Mesh.Excite.NExcite = NExcInputTot;

P.Mesh.Rotor.NDof  = NDofRotor;
P.Mesh.Stator.NDof = NDofStator;

P.Mesh.Bearing.NInput = NInputBearing;