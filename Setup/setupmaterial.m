function mat = setupmaterial(material)

if ischar(material)
    material_orig.name = material;
else
    material_orig = material;
end

mat.name = material_orig.name;
switch mat.name
    case 'steel'   
        mat.E = 200E9;
        mat.v = 0.3;
        mat.rho = 7800;
        mat.eta = 1E-4;
    case 'aluminium'
        mat.E = 70E9;
        mat.v = 0.3;
        mat.rho = 2700E9;
        mat.eta = 1E-4;
    case 'rigid'
        mat.E = Inf;
        mat.v = 0;
        mat.rho = 0;
        mat.eta = 0;
    case 'custom'
        fields = {'E','rho','v','eta'};
        for i = 1:3
            if ~isfield(material_orig,fields{i})
                error('Missing field %s from material',fields{i})
            end
        end
    otherwise
        error('Uknown material')
end


%overwrite with any userdefined properties
fields = {'E','rho','v','eta'};
for i = 1:length(fields)
    if isfield(material_orig,fields{i})
        mat.(fields{i}) = material_orig.(fields{i});
    end
end

%now compute derived properties
mat.G = mat.E/2/(1+mat.v);

mat.tensor = mat.E/(1-mat.v^2) * [ 1     mat.v       0;
                                 mat.v     1        0;
                                   0           0     0.5*(1-mat.v)];