function R = radial_offset(r,dr)
R = [1   0          dr; %w
     0   (r+dr)/r   0; %dw_dt
     0   0          1]; %dw_dr