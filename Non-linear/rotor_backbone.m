function results = rotor_backbone(P,hbm,A,O,iDofX,iDofU)
problem = rotor2problem(P);
[hbm,problem] = setuphbm(hbm,problem);
P = problem.P;

problem.res.output = 'fnl';
problem.res.iDof = iDofX;

problem.res.input = 'fe';
problem.res.iInput = iDofU;

x0 = rotor_init(hbm,problem,O,A(1));
xEnd = rotor_init(hbm,problem,O,A(end));

results = hbm_bb(hbm,problem,A(1),O,x0,A(end),O,xEnd);
results = renameStructField(results,'w','O');