function fig = plot_orbit(P,Omega,modes,omega,iNodes)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

NSpeeds = size(modes,3);
NModes = size(omega,1);
NRotors = length(P.Rotor);

col = lines(NRotors);

NRows = ceil(NModes/2);

aLocus = linspace(0,2*pi-0.3,100);

%extract the useful formation from the modeshapes
u = zeros(NRotors,NModes,NSpeeds);
v = zeros(NRotors,NModes,NSpeeds);
for k = 1:NRotors
    u(k,:,:) = modes((P.Rotor{k}.iGlobal(iNodes(k))-1)*4 + 1,:,:);
    v(k,:,:) = modes((P.Rotor{k}.iGlobal(iNodes(k))-1)*4 + 2,:,:);
end
scale = max(sqrt(abs(u).^2 + abs(v).^2),[],1);
u = u ./ repmat(scale,NRotors,1,1);
v = v ./ repmat(scale,NRotors,1,1);

fig = figure;
ax = zeros(NModes,1);
hLocus = zeros(NModes,NRotors);
hStart = zeros(NModes,NRotors);
hEnd   = zeros(NModes,NRotors);

for i = 1:NModes
    ax(i) = subplot(NRows,2,i);
    hold on
    for j = 1:NRotors
        hLocus(i,j) = plot(real(u(j,i,1) * exp(i*aLocus)), real(v(j,i,1) * exp(i*aLocus)),'color',col(j,:));
        
        hStart(i,j) = plot(ax(i), real(u(j,i,1) * exp(i*aLocus(1)))    , real(v(j,i,1) * exp(i*aLocus(1)))  ,'color',col(j,:),'marker','o');
        hEnd(i,j)   = plot(ax(i), real(u(j,i,1) * exp(i*aLocus(end)))  , real(v(j,i,1) * exp(i*aLocus(end))),'color',col(j,:),'marker','x');
    end
    
    %make the axes look pretty
    if i == 1, legend(hLocus(i,:),{P.Rotor.Name}); end
    title(ax(i),['\omega' sprintf(' = %0.2f Hz',omega(i)/(2*pi))]);
    xlabel(ax(i),'x (m)')
    ylabel(ax(i),'y (m)')
    
    %sort the axes out
    axis equal
    xlim(ax(i),[-1 1])
    ylim(ax(i),[-1 1])
end

%Rotor speed slider
hSpeedTxt = uicontrol('Style','text',...
    'Position',[335 61 60 20],...
    'String',sprintf('%d rpm',floor(Omega(1))));

hSpeedSlider = uicontrol('Style', 'slider',...
    'Min',Omega(1),'Max',Omega(end),'Value',Omega(1),...
    'Position', [207  61   120   20],'Callback',{@update_speed,fig,ax,hLocus,hStart,hEnd,hSpeedTxt,aLocus,u,v,omega,Omega}); 

uicontrol('Style','text',...
    'Position',[207 86 120 20],...
    'String','Rotor Speed');

function update_speed(obj,event,fig,ax,hLocus,hStart,hEnd,hSpeedTxt,aLocus,u,v,omega,Omega)
U = get(fig,'UserData');
[~,U.iPlot] = min(abs(Omega-obj.Value));
set(fig,'UserData',U);

NModes = size(u,2);
NRotors = size(u,1);
for i = 1:NModes
    for j = 1:NRotors
        set(hLocus(i,j),'xdata',real(u(j,i,U.iPlot) * exp(1i*aLocus))     ,'ydata',real(v(j,i,U.iPlot) * exp(1i*aLocus)));
        set(hStart(i,j),'xdata',real(u(j,i,U.iPlot) * exp(1i*aLocus(1)))  ,'ydata',real(v(j,i,U.iPlot) * exp(1i*aLocus(1))));
        set(hEnd(i,j)  ,'xdata',real(u(j,i,U.iPlot) * exp(1i*aLocus(end))),'ydata',real(v(j,i,U.iPlot) * exp(1i*aLocus(end))));
    end
    title(ax(i),['\omega' sprintf(' = %0.2f Hz',omega(i,U.iPlot)/(2*pi))]);
end

set(hSpeedTxt,'String',sprintf('%d rpm',floor(Omega(U.iPlot)*30/pi)));