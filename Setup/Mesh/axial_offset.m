function R = axial_offset(z)
R = [1   0   0   z;
     0   1  -z   0;
     0   0   1   0;
     0   0   0   1];