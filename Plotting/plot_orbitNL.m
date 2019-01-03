function [axSurf,hSurf,axOrb,hOrb] = plot_orbitNL(varargin)

if  isempty(varargin{1}) || all(ishandle(varargin{1}(:)))
    axSurf = varargin{1};
    axOrb = varargin{2};
    varargin = varargin(3:end);
    bGenPlots = isempty(axSurf);
else
    bGenPlots = 1;
end

if length(varargin) > 3
    [x,O,options,Oplot] = deal(varargin{:});
elseif length(varargin) > 2
    [x,O,options] = deal(varargin{:});
else 
    [x,O] = deal(varargin{:});
end

options = defaultmissingoptions(options);

%units for rotor speed
switch options.Ounits 
    case 'rpm'
        rad2rpm = 30/pi;
    case 'rads'
        rad2rpm = 1;
    case 'Hz'
        rad2rpm = 1/(2*pi);
end

switch options.zunits 
    case 'mm'
        m2mm = 1E3;
    case 'cm'
        m2mm = 1E2;
    case 'm'
        m2mm = 1;
    case 'mum'
        m2mm = 1E6;
        options.zunits = '\mu m';
end

NPlot = length(options.iPlot);
X = rad2rpm*repmat(O(:)',size(x,3),1,NPlot);
Y = m2mm*permute(x(options.iPlot*4-3,:,:),[3 2 1]);
Z = m2mm*permute(x(options.iPlot*4-2,:,:),[3 2 1]);

if bGenPlots
    figure
    for i = 1:NPlot
        axSurf(i) = subplot(1,NPlot,i);
        hold on
        
        da = daspect;
        da(3) = da(2);
        daspect(axSurf(i), da);
        xlabel(axSurf(i), sprintf('%s (%s)', options.Olabel, options.Ounits))
        ylabel(axSurf(i), sprintf('x (%s)', options.zunits))
        zlabel(axSurf(i), sprintf('y (%s)', options.zunits))

        set(axSurf(i),'yscale',options.zscale,'zscale',options.zscale);
        
    end
end

hSurf = zeros(NPlot,1);
for i = 1:NPlot
    hSurf(i) = surf(axSurf(i),X(:,:,i),Y(:,:,i),Z(:,:,i));
end

if length(varargin) > 3
    if bGenPlots
        figure
        for i = 1:length(Oplot)
            axOrb(i) = subplot(1,length(Oplot),i);
            hold on
            xlabel(axOrb(i), sprintf('x (%s)', options.zunits))
            ylabel(axOrb(i), sprintf('y (%s)', options.zunits))
            title(axOrb(i), sprintf('%s = %0.3f %s', options.Olabel, Oplot(i)/2/pi, options.Ounits))
            axis equal
        end
    end
    
    temp = permute(interp1(O,permute(x,[2 1 3]),Oplot),[2 3 1]);
    X = m2mm*temp(options.iPlot*4-3,:,:);
    Y = m2mm*temp(options.iPlot*4-2,:,:);
    
    hOrb = zeros(length(Oplot),1);
    for i = 1:length(Oplot)
        hOrb(i) = plot(axOrb(i),X(:,:,i)',Y(:,:,i)');
    end
end

function opt = defaultmissingoptions(opt)
f = {'Oscale','Ounits','Olabel','zunits',  'zscale'};
d = {'log'   ,'rpm'   ,'\Omega',  'm'   ,  'linear'};
for i = 1:length(f)
    if ~isfield(opt,f{i})
        opt.(f{i}) = d{i};
    end
end