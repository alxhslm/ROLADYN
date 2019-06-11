function D = dim2mass(D)

[Vol] = disc_properties(1,D.R(1),D.R(2),D.t);

D.Material.rho = D.m/Vol;

[D.m,D.Id,D.Ip]  = disc_properties(D.Material.rho,D.R(1),D.R(2),D.t);




