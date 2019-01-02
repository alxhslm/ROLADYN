function REB_analytical
uz = sym('uz');
ur = sym('ur');

alpha = sym('alpha');
K = sym('K');
A = sym('A');
Q = sym('Q');

%displacements - inner
db_duz =  sin(alpha);
db_dur =  cos(alpha);

%contact angles - outer
da_duz =  cos(alpha)/A;
da_dur = -sin(alpha)/A;

Fuz = Q*sin(alpha); %z
Fur = Q*cos(alpha); %r

Fu = [Fuz; Fur];
u = [uz; ur];

db_du = [db_duz; db_dur];
da_du = [da_duz; da_dur];

for i = 1:2
    for j = 1:2
        Juu(i,j) = diff(Fu(i),Q)*K*db_du(j) + diff(Fu(i),alpha)*da_du(j);
    end
end

% turn symbolic jacobian into a new function valJacobA.m
filename = [onedriveroot strrep('\Documents\MATLAB\ROLADYN\Bearings\REB\REB_harris_simple.m','\',filesep)];
matlabFunction(Juu,'file',filename,'vars',[K A alpha Q],'Optimize',true);

find_and_replace(filename,'\[2,2\]','[size(K),2,2]');

edit(filename);
