function results = rotor_frf_nonlin(P,hbm,A,O,x0,w)
%FRF_NONLIN Summary of this function goes here
%   Detailed explanation goes here

problem = rotor2problem(P);
problem.iDofPlot = P.Model.iDofPlot;
problem.name = P.Model.name;
[hbm,problem] = setuphbm(hbm,problem);

P = problem.P;
if nargin < 4 || isempty(x0)
    x0 = rotor_init_ode(hbm,P,O(1),A);
end

problem.X0scale = P.Model.x0scale;
problem.Xhscale = P.Model.xhscale;

problem.P.Bearing{1}.Params{2}.Options.bSaturate = 1;

if nargin>5
    %Asynchronous
    hWB = waitbar(0, 'Rotation Speed');
    for i = 1:length(O)
%         x0 = rotor_init3d(hbm,problem,[O(i),w(1)],A);
        results{i} = hbm_frf_3d(hbm,problem,O(i),w(1),w(end),A,x0); 
        results{i}.O = O(i);
    end
    close(hWB);
else
    %Synchronous so w == O
    %generate a good initial guess
%     x0 = rotor_init3d(hbm,P,O(1),A);
    results = hbm_frf(hbm,problem,O(1),O(end),A,x0);
    results.O = results.w;
    results = rmfield(results,'w');
    results.x0 = x0;
end