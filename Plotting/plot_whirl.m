function [ax,hSpeedSlider] = plot_whirl(P,Omega,modes,omega,kappa)
%PLOT_WHIRL plots the deflected shape of the rotors for each mode in 3D

NModes = size(omega,1);
NRotors = length(P.Rotor);
NBearings = length(P.Bearing);

col = lines(NRotors);

NCols = floor(sqrt(NModes));
NRows = ceil(NModes/NCols);

aLocus = linspace(0,2*pi,20);

modes = mtimesx(P.Model.A,modes);

%extract the useful formation from the modeshapes
NNodes = [];
iNodesIn = [];
for i = 1:length(P.Rotor)
    NNodes(i) = length(P.Rotor{i}.Nodes);
    iNodesIn = [iNodesIn P.Rotor{i}.iGlobal(1:(4*NNodes(i)))];
end
u = modes(iNodesIn(1:4:end),:,:);
v = modes(iNodesIn(2:4:end),:,:);
umax = max(abs(u),[],1);
vmax = max(abs(v),[],1);
% scale = 1./(repmat(max([umax;vmax],[],1),size(modes,1),1,1) + eps);
scale = 10;
modes = modes .* (scale + eps);

fig = figure('Name','Modeshapes');
ax = zeros(NModes,1);

han.Circle = {};
han.Line = {};
han.Shaft = {};
han.Disc = {};
han.Bearing = {};
han.BearingCircle = {};
for i = 1:NModes
    ax(i) = subplot(NRows,NCols,i);
%     ax(i) = tight_subplot(NRows,NCols,i,[0.01 -0.2],0,0);
    hold on
    for j = 1:NRotors
        SRotor = P.Rotor{j}.S;
        Nodes = P.Rotor{j}.Nodes;
%         plot3(Nodes,0*Nodes,0*Nodes,'color',col(j,:));

        %circles and lines at every node
        for m = 1:length(Nodes)
            zNode = Nodes(m);
            qNode = P.Rotor{j}.SNode{m} * SRotor * modes(:,i,1);
            uNode = qNode(1); vNode = qNode(2);           

            %circles at node
            han.Circle{i,j}(m) = plot3(zNode + aLocus*0,  real(uNode*exp(1i*aLocus)),   real(vNode*exp(1i*aLocus)),'color',col(j,:),'LineStyle','-');

            %line to locus at node
            han.Line{i,j}(m)   = plot3([zNode zNode], [0 real(uNode)], [0 real(vNode)] ,'color',col(j,:),'LineStyle','none');
        end

        %shaft shape
        for k = 1:length(P.Rotor{j}.Shaft)
            Shaft = P.Rotor{j}.Shaft{k};
            SShaft = Shaft.S * SRotor;
            for m = 1:(P.Rotor{j}.Shaft{k}.Mesh.Nz-1)
                zNode = P.Rotor{j}.Shaft{k}.Mesh.z(k);
                le = Shaft.Element{m}.L;
                qe = Shaft.Element{m}.S * SShaft * modes(:,i,1);
                [zShaft,uShaft,vShaft] = shaft_disp_field(qe,le,10);
                
                %shape of shaft
                han.Shaft{i,j}{k}(m)  = plot3(zNode+zShaft, real(uShaft), real(vShaft),'color',col(j,:),'LineWidth',1.5);
            end
        end
        
        %discs
        for k = 1:length(P.Rotor{j}.Disc)
            Disc = P.Rotor{j}.Disc{k};
            zDisc = Nodes(Disc.iNode);
            
            SDisc = Disc.S * SRotor;
            qd = Disc.Root.S * SDisc * modes(:,i,1);
            ud = qd(1); vd = qd(2);
            
            %circle to denote where disc is
            han.DiscRoot{i,j}(k) = plot3(zDisc, real(ud),real(vd),'color',col(j,:),'marker','o');
                        
            %shape of disc
            qDisc = Disc.Hub.S * SDisc * modes(:,i,1);
            uDisc = qDisc(1); vDisc = qDisc(2);
            for m = 1:Disc.Mesh.Nt
                for n = 1:(Disc.Mesh.Nr-1)
                    qe = Disc.Element{m,n}.S * SDisc * modes(:,i,1);
                    [xDisc,yDisc,wDisc] = disc_disp_field(qe,Disc.Mesh.r(n),Disc.Mesh.theta(m),Disc.Mesh.dr(n),Disc.Mesh.dtheta(m),10);
                    han.Disc{i,j}{k}(m,n) = surf(zDisc+real(wDisc),xDisc+real(uDisc),yDisc+real(vDisc));
                end
            end
        end
        kap = kappa(P.Rotor{j}.Disc{1}.iNode,i,1);
        if kap > 0
            dir{j} = 'FW';
        elseif kap<0
            dir{j} = 'BW';
        elseif ~any(kap)
            dir{j} = 'n/a';
        end
    end
    
    for j = 1:NBearings
        qb = P.Bearing{j}.Sb * modes(:,i,1);
        ub = qb(1); vb = qb(2);
        zb = P.Bearing{j}.z;  
        
        han.Bearing{j} = plot3(zb, real(ub),real(vb),'color',0.5*[1 1 1],'marker','x');     
        
        han.BearingCircle{j} = plot3(zb + aLocus*0,  real(ub*exp(1i*aLocus)),   real(vb*exp(1i*aLocus)),'color',0.5*[1 1 1],'LineStyle',':');  
    end
    
    %make the axes look pretty
    if i == 1
        for j = 1:NRotors
            h(j) = plot(NaN,NaN,'color',col(j,:));
            leg{j} = P.Rotor{j}.Name;
        end
%         legend(h,leg); 
    end
%     dir = dir(1);
%     zlabel(ax(i),[sprintf('Mode %d \n',i) sprintf('%s ',dir{:}) sprintf('\n%0.1f Hz',omega(i,1)/(2*pi))],'Rotation',0,'Position',[23.7383   82.2127  -59.4555]);
    title(ax(i),[sprintf('Mode %d \n',i) '\omega' sprintf(' = %0.2f Hz',omega(i,1)/(2*pi))]);
    xlabel(ax(i),'z (m)')
    ylabel(ax(i),'x (m)')
    zlabel(ax(i),'y (m)')
    
    %sort the axes out
    axis equal
    ylim(ax(i),[-1 1])
    zlim(ax(i),[-1 1])
%     set(ax(i),'YTick',[],'ZTick',[]);
%     daspect([1 2.5 2.5]);
    view(ax(i),3)
end

% update_speed(1,ax,han.Shaft,han.Circle,han.Line,han.Disc,han.DiscRoot,han.Bearing,han.BearingCircle,aLocus,P,modes,omega,kappa);
% savetikz('modes_LS.tikz',gcf,0);
% find_and_replace('modes_LS.tikz',', align=center},',', align=center, rotate=-90},')

% update_speed(length(Omega),ax,han.Shaft,han.Circle,han.Line,han.Disc,han.DiscRoot,han.Bearing,han.BearingCircle,aLocus,P,modes,omega,kappa);
% savetikz('modes_HS.tikz',gcf,0);
% find_and_replace('modes_HS.tikz',', align=center},',', align=center, rotate=-90},')

%Rotor speed slider
han.SpeedTxt = uicontrol('Style','text',...
    'Position',[335 61 60 20],...
    'String',sprintf('%d rpm',floor(Omega(1))));

han.SpeedSlider = uicontrol('Style', 'slider',...
    'Min',Omega(1),'Max',Omega(end),'Value',Omega(1),...
    'Position', [207  61   120   20],'Callback',{@slider_clbk,fig,ax,han,aLocus,P,modes,omega,kappa,Omega}); 

uicontrol('Style','text',...
    'Position',[207 86 120 20],...
    'String','Rotor Speed');

function slider_clbk(obj,event,fig,ax,han,aLocus,P,modes,omega,kappa,Omega)
[~,iPlot] = min(abs(Omega-obj.Value));
update_plot(iPlot,ax,han,aLocus,P,modes,omega,kappa);
set(han.SpeedTxt,'String',sprintf('%d rpm',floor(Omega(iPlot)*30/pi)));

function update_plot(iPlot,ax,han,aLocus,P,modes,omega,kappa)
NModes = size(modes,2);
NRotors = length(P.Rotor);
for i = 1:NModes
    for j = 1:NRotors
        SRotor = P.Rotor{j}.S;
        
        for m = 1:length(P.Rotor{j}.Nodes)
            qNode = P.Rotor{j}.SNode{m} * SRotor * modes(:,i,iPlot);
            uNode = qNode(1); vNode = qNode(2);     
            
            %circles at node
            set(han.Circle{i,j}(m), 'ydata', real(uNode*exp(1i*aLocus)), 'zdata', real(vNode*exp(1i*aLocus)));
            
            %line to locus at node
            set(han.Line{i,j}(m), 'ydata', [0 real(uNode)], 'zdata', [0 real(vNode)]);
        end
        
        for k = 1:length(P.Rotor{j}.Shaft)
            Shaft = P.Rotor{j}.Shaft{k};
            SShaft = Shaft.S * SRotor;
            for m = 1:(P.Rotor{j}.Shaft{k}.Mesh.Nz-1)
                le = Shaft.Element{m}.L;
                qe = Shaft.Element{m}.S * SShaft * modes(:,i,iPlot);
                [~,uShaft,vShaft] = shaft_disp_field(qe,le,10);
                
                %shape of shaft
                set(han.Shaft{i,j}{k}(m), 'ydata', real(uShaft), 'zdata', real(vShaft));
            end
        end
        
        for k = 1:length(P.Rotor{j}.Disc)
            Disc = P.Rotor{j}.Disc{k};
            zDisc = P.Rotor{j}.Nodes(Disc.iNode);
            
            SDisc = Disc.S * SRotor;
            qRoot = Disc.Root.S * SDisc * modes(:,i,iPlot);
            uRoot = qRoot(1); vRoot = qRoot(2);
            
            %circle to denote where disc is
            set(han.DiscRoot{i,j}(k),'ydata', real(uRoot), 'zdata', real(vRoot));
            
            %shape of disc
            qHub = Disc.Hub.S * SDisc * modes(:,i,iPlot);
            uHub = qHub(1); vHub = qHub(2);
            for m = 1:Disc.Mesh.Nt
                for n = 1:(Disc.Mesh.Nr-1)
                    qe = Disc.Se{m,n} * SDisc * modes(:,i,iPlot);
                    [xElem,yElem,wElem] = disc_disp_field(qe,Disc.Mesh.r(n),Disc.Mesh.theta(m),Disc.Mesh.dr(n),Disc.Mesh.dtheta(m),10);
                    set(han.Disc{i,j}{k}(m,n),'ydata',xElem+real(uHub),'zdata',yElem+real(vHub),'xdata',zDisc+real(wElem));
                end
            end
        end
          
        kap = kappa(P.Rotor{j}.Disc{1}.iNode,i,iPlot);
        if kap >=0
            dir{j} = 'FW';
        elseif kap<0
            dir{j} = 'BW';
        elseif ~any(kap)
            dir{j} = 'n/a';
        end
    end
    for k = 1:length(P.Bearing)
        qb = P.Bearing{k}.Sb * modes(:,i,iPlot);
        ub = qb(1); vb = qb(2);
        
        %cross to denote where bearing is
        set(han.Bearing{k},'ydata', real(ub), 'zdata', real(vb));     
        set(han.BearingCircle{k}, 'ydata', real(ub*exp(1i*aLocus)), 'zdata', real(vb*exp(1i*aLocus)));
        
    end
    title(ax(i),['\omega' sprintf(' = %0.2f Hz',omega(i,iPlot)/(2*pi))]);
%     dir = dir(1);
%     zlabel(ax(i),[sprintf('Mode %d \n',i) sprintf('%s ',dir{:}) sprintf('\n%0.1f Hz',omega(i,iPlot)/(2*pi))],'Rotation',0,'Position',[23.7383   82.2127  -59.4555]);
end

function [z,ue,ve] = shaft_disp_field(qe,le,N)
qz = linspace(0,le,N)';
C = quad_line_shapefun('coeff',le);
P = quad_line_shapefun('fun',qz);
qu = qe([1 4 5 8]);
qv = qe([2 3 6 7]); qv([2 4]) = -qv([2 4]);
ue = P*(C\qu);
ve = P*(C\qv);
z = qz;

function [xe,ye,we] = disc_disp_field(qe,r,theta,dR,beta,N)
dtheta = linspace(0,beta,N);
dr = linspace(0,dR,N);

C = quad_rect_shapefun('coeff',beta,dR);
a = C\qe;

[qh,qz] = meshgrid(dtheta,dr);
qh = qh(:); qz = qz(:);
P = quad_rect_shapefun('fun',qh,qz);
we = reshape(P*a,N,N);

theta = theta+dtheta;
r = r+dr;

theta = theta(:)';
r = r(:);

xe = r*cos(theta);
ye = r*sin(theta);