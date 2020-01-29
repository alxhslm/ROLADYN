function U = excitation_hbm(P)
U{1} = P.Mesh.Excite.uSync;
U{2} = P.Mesh.Excite.uAsync;