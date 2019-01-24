function B = setupbearing_stiffness(B,R,O,A,x0)
for i = 1:length(B)
    for j = 1:2
        qi{i}{j} = B{i}.Ri{j} * B{i}.Si{j} * x0;
        qo{i}{j} = B{i}.Ro{j} * B{i}.So{j} * x0;
    end
end

for i = 1:length(B)
    %if connected to ground, pad the iRotor vector
    if length(B{i}.iRotor) == 1
        B{i}.iRotor = [NaN B{i}.iRotor];
        B{i}.iNode  = [NaN B{i}.iNode];
    end
    
    if ~isnan(B{i}.iRotor(1))
        Oo = R{B{i}.iRotor(1)}.Speed*O;
        Ao = R{B{i}.iRotor(1)}.Speed*A;
    else
        Oo = 0;
        Ao = 0;
    end
    
    if ~isnan(B{i}.iRotor(2))
        Oi = R{B{i}.iRotor(2)}.Speed*O;
        Ai = R{B{i}.iRotor(2)}.Speed*A;
    else
        Oi = 0;
        Ai = 0;
    end
        
    B{i} = setup_each_bearing(B{i},qi{i},qo{i},Oi,Oo,Ai,Ao);
end

function B = setup_each_bearing(B,xi,xo,Oi,Oo,Ai,Ao)

States.Oi = Oi;
States.Oo = Oo;
States.Ai = Ai;
States.Ao = Ao;
wons = Ai*0 + 1;

for j = 1:2
    States.qi = xi{j}+B.ui{j}*wons;
    States.qo = xo{j}+B.uo{j}*wons;

    if strncmp(B.Model{j},'REB',3)
        States.bSolve = 1;
        
        [Forces,Channels,Stiffness] = REB_model(B.Params{j}, States);
        B.Params{j}.F0 = mean(Forces.F,2);
        B.Params{j}.K0 = mean(Stiffness.K,3);
        B.Params{j}.C0 = mean(Stiffness.C,3);

        B.Params{j}.Qi0 = mean(Channels.Qi,2);
        B.Params{j}.Qo0 = mean(Channels.Qo,2);

        B.Params{j}.Fi0 = mean(Forces.Fi,2);
        B.Params{j}.Fo0 = mean(Forces.Fo,2);
        B.Params{j}.qi0 = States.qi;
        B.Params{j}.qo0 = States.qo;
        B.Params{j}.x0 = mean(Forces.xInt,2);

        %assemble outputs
        B.Kb{j} = B.Params{j}.K0 + kron([1 -1; -1 1], B.Params{j}.Setup.KbParallel);

        B.Cb{j} = B.Params{j}.C0 + kron([1 -1; -1 1], B.Params{j}.Setup.CbParallel);

        B.Fb{j} = B.Params{j}.F0 + B.Params{j}.KPar*[States.qi; States.qo];

        B.xInt{j} = B.Params{j}.x0;

    elseif strncmp(B.Model{j},'SFD',3)
        States.bSolve = 1;
        [Forces,~,Stiffness] = SFD_model(B.Params{j}, States);
        B.Params{j}.F0 = Forces.F;
        B.Params{j}.K0 = Stiffness.K;
        B.Params{j}.C0 = Stiffness.C;
        B.Params{j}.M0 = Stiffness.M;
        B.Params{j}.qi0 = States.qi;
        B.Params{j}.qo0 = States.qo;

        B.Kb{j} = B.Params{j}.K0 + B.Params{j}.KbSquirrel;
        B.Cb{j} = B.Params{j}.C0;
        B.Fb{j} = B.Params{j}.F0 + B.Params{j}.KSq*(States.qi-States.qo);
    else
        Kb = max(-1E20,min(B.Kb{j},1E20));
        B.Fb{j} = Kb*[xi{j}+B.ui{j};xo{j}+B.uo{j}];
    end
    
    %and extract x and y components
    B.Kxx{j} = B.Kb{j}(1:2,1:2);
    B.Kxy{j} = B.Kb{j}(1:2,3:4);
    B.Kyy{j} = B.Kb{j}(3:4,3:4);
    
    B.Cxx{j} = B.Cb{j}(1:2,1:2);
    B.Cxy{j} = B.Cb{j}(1:2,3:4);
    B.Cyy{j} = B.Cb{j}(3:4,3:4);
end

%we can't have infs in the final bearing stiffness matrix as this breaks
%the FE setup code, so we set to zero now. don't worry, we'll check Kxx etc
%later to find the infs, so these values won't contribute to the final
%stiffness matrix later anyway
for j = 1:2
    ii = isinf(B.Kb{j}); B.Kb{j}(ii) = 0;
    ii = isinf(B.Cb{j}); B.Cb{j}(ii) = 0;
    %         B.Kb{j}  = max(min(B.Kb{j},1E20),-1E20);
    %         B.Cb{j}  = max(min(B.Cb{j},1E20),-1E20);
end