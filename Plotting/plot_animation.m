function ax = plot_animation(P,Omega,modes,omega,iNodes)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%various dimensions
NSpeeds = size(modes,3);
NModes = size(modes,2);
NRotors = length(P.Rotor);

NCols = floor(sqrt(NModes));
NRows = ceil(NModes/NCols);

%parameters for plotting
rShaft = (NRotors:-1:1)*0.3;
col = lines(NRotors);
aShaft = [pi linspace(0,2.5*pi,100) -0.5*pi];
aLocus = linspace(0,2*pi,100);

%extract the useful formation from the modeshapes
u = zeros(NRotors,NModes,NSpeeds);
v = zeros(NRotors,NModes,NSpeeds);
for k = 1:NRotors
    u(k,:,:) = modes((P.Rotor{k}.iGlobal(iNodes(k))-1)*4 + 1,:,:);
    v(k,:,:) = modes((P.Rotor{k}.iGlobal(iNodes(k))-1)*4 + 2,:,:);
end
scale = 0.4./max(sqrt(abs(u).^2 + abs(v).^2),[],1);
u = u .* repmat(scale,NRotors,1,1);
v = v .* repmat(scale,NRotors,1,1);

fig = figure;
ax = zeros(NModes,1);
leg = cell(NRotors,1);
hShaft = zeros(NModes,NRotors);
hLocus = zeros(NModes,NRotors);
for i = 1:NModes
    ax(i) = subplot(NRows,NCols,i);
    hold on
    for j = 1:NRotors
        %plot the outline of the shaft
        hShaft(i,j) = plot(rShaft(j)*cos(aShaft),rShaft(j)*sin(aShaft),'color',col(j,:));
        
        %plot the locus of the orbit
        hLocus(i,j) = plot(real(u(j,i,1)*exp(1i*aLocus)),real(v(j,i,1)*exp(1i*aLocus)),'color',col(j,:),'LineStyle','--');
        leg{j} = sprintf('Rotor %d',j);
    end
    
    %make the axes look pretty
    if i == 1, legend(hShaft(i,:),leg); end
    title(ax(i),['\omega' sprintf(' = %0.2f Hz',omega(i,1)/(2*pi))]);
    xlabel(ax(i),'x (m)')
    ylabel(ax(i),'y (m)')
    
    %sort the axes out
    axis equal
    xlim(ax(i),[-1 1])
    ylim(ax(i),[-1 1])
end

% Save Movie button
hSaveBtn = uicontrol('Style','pushbutton',...
    'Position',[405 11 60 20],...
    'String','Save','Callback',{@save_movie,fig});

%Playback speed slider
hPlbkTxt = uicontrol('Style','text',...
    'Position',[335 11 60 20],...
    'String','1/1000');

hPlbkSlider = uicontrol('Style', 'slider',...
    'Min',1,'Max',200,'Value',100,...
    'Position', [207  11   120   20],'Callback',{@update_plbk,hPlbkTxt}); 

uicontrol('Style','text',...
    'Position',[207 36 120 20],...
    'String','Animation Speed');

%Rotor speed slider
hSpeedTxt = uicontrol('Style','text',...
    'Position',[335 61 60 20],...
    'String',sprintf('%d rpm',floor(Omega(1))));

hSpeedSlider = uicontrol('Style', 'slider',...
    'Min',Omega(1),'Max',Omega(end),'Value',Omega(1),...
    'Position', [207  61   120   20],'Callback',{@update_speed,fig,ax,hLocus,hSpeedTxt,aLocus,u,v,omega,Omega}); 

uicontrol('Style','text',...
    'Position',[207 86 120 20],...
    'String','Rotor Speed');

U.iPlot = 1;
U.tRecord = -Inf;
U.vidObj = [];
U.tLength = 10;
U.tCurrent = 0;
set(fig,'UserData',U);

tmr = timer('TimerFcn',{@update_plots,fig,ax,hShaft,hPlbkSlider,rShaft,aShaft,u,v,omega,Omega,getrotorspeeds(P.Rotor)},'Period',0.05,'ExecutionMode','fixedrate');
start(tmr)

set(fig,'CloseRequestFcn',@(obj,event)close_fig(obj,event,tmr))

function save_movie(obj,event,fig)
U = get(fig,'UserData');
answer = inputdlg('Video length (s)','Specify Video Length');
U.tLength = str2double(answer);

[file,path] = uiputfile('*.avi','Movie file name');
U.vidObj = VideoWriter([path file]);
U.vidObj.FrameRate = 20;
open(U.vidObj);
U.tRecord = U.tCurrent;
set(fig,'UserData',U);

function update_plots(obj,event,fig,ax,hShaft,hPlbkSlider,rShaft,aShaft,u,v,omega,Omega,speedratio)
U = get(fig,'UserData');
U.tCurrent = U.tCurrent + obj.Period/hPlbkSlider.Value;

NModes = size(u,2);
NRotors = size(u,1);

for i = 1:NModes
    for j = 1:NRotors
        p = exp(1i * omega(i) * U.tCurrent);
        x = real(u(j,i,U.iPlot) * p);
        y = real(v(j,i,U.iPlot) * p);
        th = aShaft + Omega(U.iPlot)*speedratio(j)*U.tCurrent;
        set(hShaft(i,j),'xdata',rShaft(j)*cos(th) + x,'ydata',rShaft(j)*sin(th) + y);
    end
end

if  U.tCurrent < (U.tRecord + U.tLength/hPlbkSlider.Value)
    writeVideo(U.vidObj,getframe(fig));
elseif ~isempty(U.vidObj)
	close(U.vidObj);
    U.vidObj = [];    
end

set(fig,'UserData',U);

drawnow

function update_speed(obj,event,fig,ax,hLocus,hSpeedTxt,aLocus,u,v,omega,Omega)
U = get(fig,'UserData');
[~,U.iPlot] = min(abs(Omega-obj.Value));
set(fig,'UserData',U);

NModes = size(u,2);
NRotors = size(u,1);
for i = 1:NModes
    for j = 1:NRotors
        set(hLocus(i,j),'xdata',real(u(j,i,U.iPlot) * exp(1i*aLocus)),'ydata',real(v(j,i,U.iPlot) * exp(1i*aLocus)));
    end
    title(ax(i),['\omega' sprintf(' = %0.2f Hz',omega(i,U.iPlot)/(2*pi))]);
end

set(hSpeedTxt,'String',sprintf('%d rpm',floor(Omega(U.iPlot)*30/pi)));

function update_plbk(obj,event,hPlbkTxt)
set(hPlbkTxt,'String',sprintf('1/%d',floor(obj.Value)));

function close_fig(fig,event,tmr)
stop(tmr);
delete(fig);