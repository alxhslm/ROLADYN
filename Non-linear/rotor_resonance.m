function results = rotor_resonance(P,hbm,A,O,iDofX,iDofU,x0)
problem = rotor2problem(P);
[hbm,problem] = setuphbm(hbm,problem);
P = problem.P;
   
problem.res.output = 'fnl';
problem.res.iDof = iDofX;

problem.res.input = 'fe';
problem.res.iInput = iDofU;

if nargin < 7 || ~isempty(x0)
    x0 = rotor_init(hbm,problem,O,A);
end

results = hbm_res(hbm,problem,O,A,x0);
results = renameStructField(results,'w','O');
