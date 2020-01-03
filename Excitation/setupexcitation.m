function P = setupexcitation(P)
E = P.Excite;
Speed = getrotorspeeds(P.Rotor);

iCount = 0;
for i = 1:length(E)
    switch E{i}.Type
        case 'unbalance'
            E{i} = setupunbalance(E{i},P.Rotor);
        case 'skew'
            E{i} = setupskew(E{i},P.Rotor);
        case 'shaker'
            E{i} = setupshaker(E{i});
        otherwise
            error('Unrecognised excitation type')
    end
    
    if ~isfield(E{i},'Name')
        E{i}.Name = E{i}.Type;
    end
    
    fields = {'K','C','M'};
    for j = 1:length(fields)
        if ~isfield(E{i},fields{j})
            E{i}.(fields{j}) = zeros(E{i}.NInput);
        end
    end
    
    E{i}.iExcite = iCount + (1:E{i}.NInput);
    iCount = iCount + E{i}.NInput;
end

P.Excite = E;

function E = setupunbalance(E,R)
essential_fields = {'iRotor', 'iDisc'};
for j = 1:length(essential_fields)
    if ~isfield(E,essential_fields{j})
        error(['Missing field ' essential_fields{j} ' in P.Excite'])
    end
end

defaultable_fields = {'m','r', 'Angle'};
val = [R{E.iRotor}.Disc{E.iDisc}.Inertia.m 0 0];
for j = 1:length(defaultable_fields)
    if ~isfield(E,defaultable_fields{j})
        E.(defaultable_fields{j}) = val(j);
    end
end

E.NInput = 2;
E.M = eye(2);
E.u = E.m * E.r.*exp(1i*(E.Angle + [0; -pi/2*sign(R{E.iRotor}.Speed)]));
if ~isfield(E,'Mode')
    E.Mode = 'sync';
end

function E = setupskew(E,R)
essential_fields = {'iRotor', 'iDisc'};
for j = 1:length(essential_fields)
    if ~isfield(E,essential_fields{j})
        error(['Missing field ' essential_fields{j} ' in P.Excite'])
    end
end

defaultable_fields = {'Skew','Angle'};
for j = 1:length(defaultable_fields)
    if ~isfield(E,defaultable_fields{j})
        E.(defaultable_fields{j}) = 0;
    end
end

E.NInput = 2;
E.M = (R{E.iRotor}.Disc{E.iDisc}.Inertia.Id - R{E.iRotor}.Disc{E.iDisc}.Inertia.Ip)*eye(2);
E.u = E.Skew.*exp(1i*(E.Angle + pi/2 + [0; -pi/2*sign(R{E.iRotor}.Speed)]));

if ~isfield(E,'Mode')
    E.Mode = 'sync';
end

function E = setupshaker(E)
essential_fields = {'iStator','Mode'};
for j = 1:length(essential_fields)
    if ~isfield(E,essential_fields{j})
        error(['Missing field "' essential_fields{j} '" in P.Excite'])
    end
end

defaultable_fields = {'Amplitude', 'Phase'};
for j = 1:length(defaultable_fields)
    if ~isfield(E,defaultable_fields{j})
        E.(defaultable_fields{j}) = zeros(2,1);
    end
end

fields = {'K','C','M'};
if ~any(isfield(E,fields))
    error('Need either "K","C" or "M"" in P.Excite for shake excitation')
end
for i = 1:length(fields)
    if ~isfield(E,fields{i})
        E.(fields{i}) = zeros(2);
    end
end

E.NInput = 2;
E.u = E.Amplitude.*exp(1i*E.Phase);