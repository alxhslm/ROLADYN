function plot_REB_animation(REB,varargin)

thisArg = getargs(varargin,1);
fig = plot_REB(REB,thisArg{:});

U = get(fig,'UserData');
U.iTimer = 1;
U.incTimer = 1;
set(fig,'UserData',U);

tmr = timer('TimerFcn',{@update_plots,fig,REB,varargin},'Period',0.05,'ExecutionMode','fixedrate');
set(fig,'CloseRequestFcn',@(obj,event)close_fig(obj,event,tmr))

start(tmr)

function update_plots(obj,event,fig,REB,args)
U = get(fig,'UserData');

thisArg = getargs(args,U.iTimer);
plot_REB(fig,REB,thisArg{:});
drawnow
if U.iTimer == size(args{1},2)
    U.incTimer = -1;
elseif U.iTimer == 1
    U.incTimer = 1;
end
U.iTimer = U.iTimer + U.incTimer;

set(fig,'UserData',U);

function arg = getargs(args,i)
for j = 1:length(args)
    arg{j} = args{j}(:,i);
end

function close_fig(fig,event,tmr)
stop(tmr);
delete(fig);