function results = rotor_solve(P,hbm,A,O,x0,w)
problem = rotor2problem(P);
[hbm,problem] = setuphbm(hbm,problem);
P = problem.P;

%find the base frequencies
if nargin>5
    %Asynchronous
    hbm.harm.rFreqRatio = [O w]./O;
    w0 = [O w];
else
    %Synchronous so w == O
    hbm.harm.rFreqRatio = [1 0];
    w0 = O;
end

%generate a good initial guess
if isempty(x0)
    x0 = rotor_init3d(hbm,problem,w0,A);
end

results = hbm_solve(hbm,problem,w0,A,x0);   
results = renameStructField(results,'w','O');