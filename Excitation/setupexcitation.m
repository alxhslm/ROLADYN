function P = setupexcitation(P)
E = P.Excite;
Speed = getrotorspeeds(P.Rotor);

iCount = 0;
for i = 1:length(E)
    if strcmp(E{i}.Name,'unbalance')
        E{i} = setupunbalance(E{i},Speed);
        E{i}.M = P.Rotor{E{i}.iRotor}.Disc{E{i}.iDisc}.M;
        E{i}.C = zeros(4);
        E{i}.K = zeros(4);
        E{i}.NInput = 4;
    elseif  strcmp(E{i}.Name,'ground')
        E{i} = setupground(E{i});
        E{i}.M = zeros(4);
        E{i}.K = P.Bearing{E{i}.iBearing}.Kb;
        E{i}.C = P.Bearing{E{i}.iBearing}.Cb;
        E{i}.NInput = 8;
    else
        E{i}.NInput = 0;
    end
    
    E{i}.iExcite = iCount + (1:E{i}.NInput);
    iCount = iCount + E{i}.NInput;
    E{i}.fun = str2func(['exc_', E{i}.Name]);
end

P.Excite = E;

function E = setupunbalance(E,Speed)
essential_fields = {'iRotor', 'iDisc'};
for j = 1:length(essential_fields)
    if ~isfield(E,essential_fields{j})
        error(['Missing field ' essential_fields{j} ' in P.Excite'])
    end
end

defaultable_fields = {'eUnbalance','aSkew', 'aUnbalancePh','aSkewPh'};
for j = 1:length(defaultable_fields)
    if ~isfield(E,defaultable_fields{j})
        E.(defaultable_fields{j}) = 0;
    end
end

if length(E.aUnbalancePh) == 1
    E.aUnbalancePh = E.aUnbalancePh + [0;-pi/2*sign(Speed(E.iRotor))];
end

if length(E.eUnbalance) == 1
    E.eUnbalance = repmat(E.eUnbalance,2,1);
end

if length(E.aSkewPh) == 1
    E.aSkewPh = E.aSkewPh + [0;-pi/2*sign(Speed(E.iRotor))];
end

if length(E.aSkew) == 1
    E.aSkew = repmat(E.aSkew,2,1);
end

E.u = [E.eUnbalance.*exp(1i*E.aUnbalancePh);
       E.aSkew     .*exp(1i*E.aSkewPh)];

function E = setupground(E)
essential_fields = {'iBearing'};
for j = 1:length(essential_fields)
    if ~isfield(E,essential_fields{j})
        error(['Missing field ' essential_fields{j} ' in P.Excite'])
    end
end

defaultable_fields = {'uGnd', 'aGndPh'};
for j = 1:length(defaultable_fields)
    if ~isfield(E,defaultable_fields{j})
        E.(defaultable_fields{j}) = zeros(4,1);
    end
end

for j = 1:length(defaultable_fields)
    if ~iscell(E.(defaultable_fields{j}))
        E.(defaultable_fields{j}) = {E.(defaultable_fields{j}), zeros(4,1)};
    end
end

E.u = [];
for j = 1:2
    E.u = [E.u; E.uGnd{j}.*exp(1i*E.aGndPh{j})];
end