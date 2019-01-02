function varargout = trig_rect_shapefun(output,h,z)
zer = 0*h;
won = zer+1;

if ~iscell(output)
    output = {output};
end

for i = 1:length(output)
    switch  output{i}
        case 'coeff'
        varargout{i} = [trig_rect_shapefun('fun',0,0);
                        trig_rect_shapefun('jac',0,0);
                        trig_rect_shapefun('fun',0,z);
                        trig_rect_shapefun('jac',0,z);
                        trig_rect_shapefun('fun',h,z);
                        trig_rect_shapefun('jac',h,z);
                        trig_rect_shapefun('fun',h,0);
                        trig_rect_shapefun('jac',h,0);];
        case 'fun'
            varargout{i} =  [won z z.^2 z.^3 sin(h) z.*sin(h) z.^2*sin(h) z.^3*sin(h) cos(h) z.*cos(h) z.^2*cos(h) z.^3*cos(h)];
        case 'jac'   
            varargout{i} =   [zer   zer   zer    zer    cos(h)   z.*cos(h)   z.^2*cos(h)   z.^3*cos(h) -sin(h) -z.*sin(h)  -z.^2*sin(h)    z.^3*sin(h);
                              zer   won   2*z   3*z.^2   zer      sin(h)     2*z*sin(h)   3*z.^2*sin(h)   zer    cos(h)      2*z*cos(h)   3*z.^2*cos(h)];
        case 'hess' 
            varargout{i}  = [zer   zer   zer    zer   -sin(h)  -z.*sin(h)  -z.^2*sin(h)  -z.^3*sin(h) -cos(h) -z.*cos(h)  -z.^2*cos(h)    z.^3*cos(h);
                             zer   zer    2     6*z     zer       zer        2*sin(h)      6*z*sin(h)   zer      zer         2*cos(h)     6*z*cos(h);
                             zer   zer   zer    zer     zer      cos(h)    2*z*cos(h)   3*z.^2*cos(h)   zer    -sin(h)     -2*z*sin(h)  3*z.^2*sin(h)];
    end
end