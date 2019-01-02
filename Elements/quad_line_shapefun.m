function varargout = quad_line_shapefun(output,z)
zer = 0*z;
won = zer+1;

if ~iscell(output)
    output = {output};
end

for i = 1:length(output)
    switch  output{i}
        case 'coeff'
            varargout{i} = [quad_line_shapefun('fun',0);
                            quad_line_shapefun('jac',0)
                            quad_line_shapefun('fun',z);
                            quad_line_shapefun('jac',z)];
        case 'fun'
            varargout{i} =  [won    z   z.^2    z.^3];
        case 'jac'   
            varargout{i} =  [zer   won    z     3*z.^2];
        case 'hess' 
            varargout{i} =  [zer   zer   won    6*z];
    end
end