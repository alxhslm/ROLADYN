function varargout = rotor_solve(P,hbm,A,O,x0,w)
%ROTOR_SOLVE Summary of this function goes here
%   Detailed explanation goes here

problem = rotor2problem(P);
[hbm,problem] = setuphbm(hbm,problem);
P = problem.P;

if isempty(x0)
    if nargin>5
        x0 = rotor_init3d(hbm,problem,[O,w],A);
    else
        x0 = rotor_init3d(hbm,problem,O,A);
    end
end

if nargin>5
    %Asynchronous
    [varargout{1:nargout}] = hbm_solve_3d(hbm,problem,O,w,A,x0);   
else
    %Synchronous so w == O
    %generate a good initial guess
    [varargout{1:nargout}] = hbm_solve(hbm,problem,O,A,x0);   
end