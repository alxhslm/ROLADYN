function [ax_o,ax_z,hLine_o,hLine_z] = plot_campbell(varargin) 
if  isempty(varargin{1}) || all(ishandle(varargin{1}))
    ax_o = varargin{1};
    ax_z = varargin{2};
    varargin = varargin(3:end);
    bGenPlots = isempty(ax_o);
else
    bGenPlots = 1;
end

[P,Omega,omega,zeta,kappa,iNodes,options] = deal(varargin{:});

options = defaultmissingoptions(options);

NModes = size(omega,1);
NRows = ceil(NModes/2);

if bGenPlots
    fo = figure('Name','Campbell Diagram: Natural Frequenies');
    if options.bSeperateModes
        for i = 1:NModes
            ax_o(i) = subplot(NRows,2,i);
            title(sprintf('Mode %d',i));
            hold on
        end
    else
        ax_o = subplot(1,1,1);
        hold on
        ax_o(1:NModes) = ax_o;
    end

    fz = figure('Name','Campbell Diagram: Damping Coefficients');
    if options.bSeperateModes
        for i = 1:NModes
            ax_z(i) = subplot(NRows,2,i);
            title(sprintf('Mode %d',i));
            hold on
        end
    else
        ax_z = subplot(1,1,1);
        hold on
        ax_z(1:NModes) = ax_z;
    end
else
    if numel(ax_o) == 1
         ax_o(1:NModes) = ax_o;
    end
    if numel(ax_z) == 1
         ax_z(1:NModes) = ax_z;
    end
end

col = lines(NModes);

leg_modes = cell(NModes,1);
hLine_o   = zeros(NModes,1);
hLine_z   = zeros(NModes,1);

if options.bAsymmetric
    freq = 0*repmat(omega,2,1);
    for j = 1:NModes
        freq(2*j-1,:) = abs(Omega + omega(j,:));
        freq(2*j,:) = abs(Omega - omega(j,:));
    end
else
    freq = omega;
end

switch options.xunits 
    case 'rpm'
        rad2rpm = 30/pi;
    case 'rads'
        rad2rpm = 1;
    case 'Hz'
        rad2rpm = 1/(2*pi);
end

switch options.yunits{1}
    case 'rpm'
        rad2hz = 30/pi;
    case 'rads'
        rad2hz = 1;
    case 'Hz'
        rad2hz = 1/(2*pi);
end

switch options.yunits{2}
    case ''
        frac2percent = 1;
    case '%'
        frac2percent = 100;
end


NRotor = length(P.Rotor);
dir = cell(NRotor,1);
for j = 1:NModes
    if options.bAsymmetric
        hLine_o(j) = plot(ax_o(j),Omega*rad2rpm,freq(2*j-1,:)*rad2hz,'color',col(j,:));
        plot(ax_o(j),Omega*rad2rpm,freq(2*j,:)*rad2hz,'color',col(j,:));
    else
        hLine_o(j) =  plot(ax_o(j),Omega*rad2rpm,freq(j,:)*rad2hz,'color',col(j,:));
    end
    hLine_z(j) = plot(ax_z(j),Omega*rad2rpm,zeta(j,:)*frac2percent,'color',col(j,:));

    for i = 1:NRotor
        kap = kappa(P.Rotor{i}.iGlobal(iNodes(i)),j,2:end);
        if mostly(kap>0,0.95)
            dir{i} = 'FW';
        elseif mostly(kap<0,0.95)
            dir{i} = 'BW';
        elseif ~any(kap)
            dir{i} = 'n/a';
        else
            %the direction changes
            if kap(1) > 0
                dir{i} = 'FW\rightarrowBW';
            else
                dir{i} = 'BW\rightarrowFW';
            end
        end
    end
%      leg_modes{j} = [sprintf('Mode %d ',j) '(' sprintf('%s;',dir{:}) ')'];
     leg_modes{j} = sprintf('Mode %d (%s)',j,dir{1});
end

if bGenPlots
    hLeg = hLine_o;
    leg = leg_modes;
    style = {'-','-','-.',':'};
    
    for j = 1:NModes
        if j == 1
            for k = 1:length(options.engineorder)
                h = plot(ax_o(j),Omega*rad2rpm,Omega*options.engineorder(k)*rad2hz,'LineStyle',style{k},'color','k');
                if ~isempty(options.engineorderlabel{k})
                    hLeg(end+1) = h;
                    leg{end+1} = options.engineorderlabel{k};
                end
            end
        end

        xlim(ax_o(j),(mean(Omega) + max(range(Omega),1)*[-0.5 0.5])*rad2rpm);
        ylim(ax_o(j),[0 max(freq(:))]*rad2hz);
        xlabel(ax_o(j),getlabel(options,'x'));
        ylabel(ax_o(j),getlabel(options,'y1'))

        xlim(ax_z(j),(mean(Omega) + max(range(Omega),1)*[-0.5 0.5])*rad2rpm);
        ylim(ax_z(j),(mean(zeta(:)) + max(range(zeta(:)),1E-3)*[-0.75 0.75])*frac2percent);
        xlabel(ax_z(j),getlabel(options,'x'));
        ylabel(ax_z(j),getlabel(options,'y2'));
        
    end

    if ~options.bSeperateModes
        ax_o = ax_o(1);
        ax_z = ax_z(1);
        legend(ax_o,hLeg,leg,'Location','EastOutside');
        legend(ax_z,hLine_z,leg_modes,'Location','EastOutside');
        linkaxes([ax_o ax_z],'x');      
    else
        linkaxes([ax_o(:); ax_z(:)],'x');    
        linkaxes(ax_o(:),'y'); 
        linkaxes(ax_z(:),'y'); 
    end

end

function opt = defaultmissingoptions(opt)
f = {'bSeperateModes','bAsymmetric','bUnitsInLabels','engineorder','xunits','yunits','xlabel','ylabel'};
d = {false,false,true,1:4, 'rpm',{'Hz',''},'\Omega',{'\omega','\zeta'}};
for i = 1:length(f)
    if ~isfield(opt,f{i})
        opt.(f{i}) = d{i};
    end
end
if ~isfield(opt,'engineorderlabel')
    for k = 1:length(opt.engineorder)
        opt.engineorderlabel{k} = sprintf('%d\\Omega',opt.engineorder(k));
    end
end

function b = mostly(x,tol)
%have some tolerance to allow for small bugs in the continuation
b = sum(x)/length(x) > tol;

function s = getlabel(options,t)
switch t
    case 'x'
            s = options.xlabel;
            u = options.xunits;
    case 'y1'
            s = options.ylabel{1};
            u = options.yunits{1};
    case 'y2'
            s = options.ylabel{2};
            u = options.yunits{2};
end

if options.bUnitsInLabels
    s = sprintf('%s (%s)', s, u);
end