function init = rotor_init_ode(harm,P,O,A)
%work out frequencies
w0 = O*harm.rFreqRatio;

T = 2*pi/O;
N = 100;
dt = T/N;
time = 0:dt:50*T;

[t0, x0] = rotor_ode(P, time, O, w0(2), A);

iKeep = find(t0 < floor(t0(end)/T)*T,1,'last');
Nkeep = floor(0.5*t0(end)/T);
iKeep = iKeep + 1 + (-Nkeep*N:-1);
t = t0(iKeep);
x = x0(iKeep,:);
X = REB_freq_comp(t,x,harm,P,w0);

init = struct('X',X,'x',x,'t',t,'O',O);