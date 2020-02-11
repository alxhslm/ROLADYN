function R = axial_offset(z)
if isnan(z)
    R = eye(4);
    return;
end
R = [1   0   0   z;
     0   1  -z   0;
     0   0   1   0;
     0   0   0   1];