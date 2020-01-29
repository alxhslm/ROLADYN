function [ax_mag,ax_ph,h_mag,h_ph] = plot_frf(varargin)
%PLOT_FRF Function to plot FRFs
%   Usage: 
%       [ax_mag,ax_ph,h_mag,h_ph] = plot_frf([ax_mag],[ax_ph],P,Q,O,[w],options)

if  isempty(varargin{1}) || all(ishandle(varargin{1}(:)))
    ax_mag = varargin{1};
    ax_ph = varargin{2};
    varargin = varargin(3:end);
    bGenPlots = isempty(ax_mag);
else
    bGenPlots = 1;
end

if length(varargin) > 4
    [P,Q,O,w,options] = deal(varargin{:});
    bAsync = true;
else
    [P,Q,O,options] = deal(varargin{:});
    bAsync = false;
end

if ~iscell(Q)
    Q = {Q};
end

NDof = size(Q{1},1);
NInput = length(Q);
options = defaultmissingoptions(options,NDof,NInput);

NOutput = length(options.iPlot);

if ~iscell(O)
    O = repmat({O},NInput,1);
end

output = index2dofnames(P,options.iPlot);

%units of asynchronous frequency
if bAsync
    switch options.wunits 
        case 'rpm'
            rad2Hz = 30/pi;
        case 'rads'
            rad2Hz = 1;
        case 'Hz'
            rad2Hz = 1/(2*pi);
    end
end

%units for rotor speed
switch options.Ounits 
    case 'rpm'
        rad2rpm = 30/pi;
    case 'rads'
        rad2rpm = 1;
    case 'Hz'
        rad2rpm = 1/(2*pi);
end

switch options.zunits{2} 
    case 'rad'
        rad2deg = 1;
    case 'deg'
        rad2deg = 180/pi;
end


%Initialise some figures and axes
if bGenPlots
    mag = figure('Name','FRF: Magnitude');
    ax_mag = gobjects(NOutput,NInput);
    for i = 1:NOutput
        for j = 1:NInput
            ax_mag(i,j) = subplot(NOutput,NInput,(i-1)*NInput + j);
            hold on
        end
    end

    pha = figure('Name','FRF: Phase');
    ax_ph = gobjects(NOutput,NInput);
    for i = 1:NOutput
        for j = 1:NInput
            ax_ph(i,j) = subplot(NOutput,NInput,(i-1)*NInput + j);
            hold on
        end
    end
end

h_mag = zeros(NOutput,NInput);
h_ph  = zeros(NOutput,NInput);

for j = 1:NInput
    q = mtimesx(P.Model.A,Q{j});   
    
    for i = 1:NOutput
        if bAsync %async - 3d surfaces
            m = abs(squeeze(q(options.iPlot(i),:,:))');
            p = phase(squeeze(q(options.iPlot(i),:,:))');
            if strcmpi(options.zunits{1},'dB')
                m = 20*log10(m);
            end
            
            h_mag(i,j) = surf(ax_mag(i,j),w*rad2Hz,O{j}*rad2rpm,  m,'EdgeColor','none','FaceColor','Interp');
            h_ph(i,j) = surf(ax_ph(i,j) ,w*rad2Hz,O{j}*rad2rpm,rad2deg*p,'EdgeColor','none','FaceColor','Interp');
        else %sync - just plots
            m = abs(q(options.iPlot(i),:));
            p = phase(q(options.iPlot(i),:));
            if strcmpi(options.zunits{1},'dB')
                m = 20*log10(m);
            end
            h_mag(i,j) = plot(ax_mag(i,j),O{j}*rad2rpm,m);
            h_ph(i,j) = plot(ax_ph(i,j) ,O{j}*rad2rpm,rad2deg*p);
        end

        if bGenPlots               
            if bAsync %async - use z label
                set(ax_mag(i,j),'XScale',options.wscale);
                set(ax_ph(i,j) ,'XScale',options.wscale);
                set(ax_mag(i,j),'YScale',options.Oscale);
                set(ax_ph(i,j) ,'YScale',options.Oscale);
                set(ax_mag(i,j),'ZScale',options.zscale);
                if i == NOutput
                    xlabel(ax_mag(i,j),sprintf('%s (%s)', options.wlabel, options.wunits))
                    xlabel(ax_ph(i,j) ,sprintf('%s (%s)', options.wlabel, options.wunits))
                end
                ylabel(ax_mag(i,j),sprintf('%s (%s)', options.Olabel, options.Ounits))
                ylabel(ax_ph(i,j) ,sprintf('%s (%s)', options.Olabel, options.Ounits))
                if j == 1
                    zlabel(ax_mag(i,j),sprintf('|%s| (%s)', output{i}, options.zunits{1}))
                    zlabel(ax_ph(i,j), sprintf('\\angle %s (%s)', output{i}, options.zunits{2}))
                end
            else %sync - use the y label
                set(ax_mag(i,j),'XScale',options.Oscale);
                set(ax_ph(i,j) ,'XScale',options.Oscale);
                set(ax_mag(i,j),'YScale',options.zscale);
                if i == NOutput
                    xlabel(ax_mag(i,j),sprintf('%s (%s)', options.Olabel, options.Ounits))
                    xlabel(ax_ph(i,j) ,sprintf('%s (%s)', options.Olabel, options.Ounits))
                end
                if j == 1
                    ylabel(ax_mag(i,j),sprintf('|%s| (%s)', output{i}, options.zunits{1}))
                    ylabel(ax_ph(i,j),sprintf('\\angle %s (%s)', output{i}, options.zunits{2}))
                end
            end

            if i == 1
                title(ax_mag(i,j),options.input{j})
                title(ax_ph(i,j),options.input{j})
            end
        end
    end
end

if bGenPlots
    if bAsync %async - link the camera views
        hLink = linkprop([ax_mag(:);ax_ph(:)],'View');
        set(mag,'UserData',hLink);
    else %sync - match the xlim
        linkaxes([ax_mag(:);ax_ph(:)],'x');
    end
end

function opt = defaultmissingoptions(opt,NDof,NInput)
f = {'wscale','wunits','wlabel','Oscale','Ounits','Olabel','zunits',   'zscale','iPlot','input'};
d = {'log'   , 'Hz'   ,'\omega','log'   ,'rpm'   ,'\Omega',{'dB','deg'},'linear',1:NDof,cellsprintf('Input %d',num2cell(1:NInput))};
for i = 1:length(f)
    if ~isfield(opt,f{i})
        opt.(f{i}) = d{i};
    end
end