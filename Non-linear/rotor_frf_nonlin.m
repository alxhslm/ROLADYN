function results = rotor_frf_nonlin(P,hbm,A,O,x0,xEnd)
problem = rotor2problem(P);
[hbm,problem] = setuphbm(hbm,problem);
P = problem.P;

%generate a good initial guess
if nargin < 4 || isempty(x0)
    x0 = rotor_init3d(hbm,P,O(1),A);
end
if nargin < 5 || isempty(xEnd)
    xEnd = rotor_init3d(hbm,P,O(end),A);
end

%now compute the frf
results = hbm_frf(hbm,problem,A,O(1),x0,O(end),xEnd);
results = renameStructField(results,'w','O');
