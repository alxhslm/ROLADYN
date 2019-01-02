function varargout = quad_rect_shapefun(output,h,z)
zer = 0*h;
won = zer+1;

if ~iscell(output)
    output = {output};
end

for i = 1:length(output)
    switch  output{i}
        case 'coeff'
        varargout{i} = [quad_rect_shapefun('fun',0,0);
                        quad_rect_shapefun('jac',0,0);
                        quad_rect_shapefun('fun',0,z);
                        quad_rect_shapefun('jac',0,z);
                        quad_rect_shapefun('fun',h,z);
                        quad_rect_shapefun('jac',h,z);
                        quad_rect_shapefun('fun',h,0);
                        quad_rect_shapefun('jac',h,0);];
        case 'fun'
            varargout{i} =  [won    h      z     h.^2    h.*z    z.^2     h.^3     (h.^2).*z     h.*(z.^2)    z.^3   (h.^3).*z    h.*(z.^3)];
        case 'jac'   
            varargout{i} = [zer   won    zer    2*h       z      zer     3*(h.^2)   2*h.*z       z.^2        zer    3*(h.^2).*z    z.^3;
                            zer    zer   won    zer       h      2*z     zer        h.^2        2*h.*z      3*z.^2    h.^3      3*h.*z.^2];
        case 'hess' 
            varargout{i}  = [zer    zer    zer    2*won    zer      zer       6*h      2*z       zer      zer      6*h.*z        zer;   % d2w/dh2
                             zer    zer    zer     zer     zer     2*won      zer      zer       2*h      6*z       zer         6*h.*z; % d2w/dz2
                             zer    zer    zer     zer     won     zer        zer      2*h       2*z      zer      3*h.^2       3*z.^2];  % d2w/dzdh
    end
end