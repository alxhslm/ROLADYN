function varargout = plot_system(P,bLabel)

if nargin < 2
    bLabel = 1;
end

col = lines(length(P.Rotor));

figure;
hold on
h.Axis = gca;
h.Figure = gcf;

rNodes = {};

for i = 1:length(P.Rotor)  
    %plot the shafts
    rNodes{i} = 0*P.Rotor{i}.Nodes;
    
    for j = 1:length(P.Rotor{i}.Shaft)
        for k = 1:length(P.Rotor{i}.Shaft{j}.le)
            [X,Y,Z] = plot_cylinder(P.Rotor{i}.Shaft{j}.ri, P.Rotor{i}.Shaft{j}.ro, P.Rotor{i}.Shaft{j}.ze(k), P.Rotor{i}.Shaft{j}.le(k));
            h.Rotor{i}.Shaft{j} = surf(Z,X,Y,'edgecolor','k','facecolor',col(i,:),'facealpha',0.3);
            rNodes{i}(P.Rotor{i}.Shaft{j}.iNodes(k)) = max(P.Rotor{i}.Shaft{j}.ro,rNodes{i}(P.Rotor{i}.Shaft{j}.iNodes(k)));
        end
    end
    
    %and the discs
    for j = 1:length(P.Rotor{i}.Disc)
        for k = 1:P.Rotor{i}.Disc{j}.NSegments
            [X,Y,Z] = plot_cylinder(P.Rotor{i}.Disc{j}.R(k),P.Rotor{i}.Disc{j}.R(k+1),P.Rotor{i}.Disc{j}.z, P.Rotor{i}.Disc{j}.t(k));
            h.Rotor{i}.Disc{j}(k) = surf(Z,X,Y,'edgecolor','k','facecolor',col(i,:));
        end
        rNodes{i}(P.Rotor{i}.Disc{j}.iNode) = max(P.Rotor{i}.Disc{j}.R(end),rNodes{i}(P.Rotor{i}.Disc{j}.iNode));
    end 
    
    Rotor_Names{i} = P.Rotor{i}.Name;
end

%now plot the bearings
for i = 1:length(P.Bearing)
    roBearing = 0;
    riBearing = Inf;
    for j = 1:2
        iRotor = P.Bearing{i}.iRotor(j);
        if ~isnan(iRotor)
            for k = 1:length(P.Rotor{iRotor}.Shaft)
                if any(P.Rotor{iRotor}.Shaft{k}.iNodes ==  P.Bearing{i}.iNode(j))
                    if j == 1 %then inner rotor
                        riBearing = min(riBearing,P.Rotor{iRotor}.Shaft{k}.ro);
                    else %otherwise this is the outer rotor
                        roBearing = max(roBearing,P.Rotor{iRotor}.Shaft{k}.ri);
                    end
                end
            end
        end
    end
    
    %if the bearing is connected to ground, then arbitrarily make it 0.01m
    % thick
    if isnan(P.Bearing{i}.iRotor(2))
        roBearing = riBearing + 0.01;
        rNodes{P.Bearing{i}.iRotor(1)}(P.Bearing{i}.iNode(1)) = max(roBearing,rNodes{P.Bearing{i}.iRotor(1)}(P.Bearing{i}.iNode(1)));
    end
    
    [X,Y,Z] = plot_cylinder(riBearing,roBearing,P.Bearing{i}.z,P.Bearing{i}.t);
    h.Bearing{i}.Obj = surf(Z,X,Y,'edgecolor','k','facecolor',[0.5 0.5 0.5]);
    
end

rText = max(cellfun(@max,rNodes)) + 0.1;

if bLabel
    %now plot the labels for the bearings
    for i = 1:length(P.Bearing)
        rPlot = 0;
        for j = 1:2
            if ~isnan(P.Bearing{i}.iRotor(j))
                rPlot = max(rPlot,rNodes{P.Bearing{i}.iRotor(j)}(P.Bearing{i}.iNode(j)));
            end
        end

        h.Bearing{i}.Label = text(P.Bearing{i}.z,-rText,0,sprintf('Bearing %d',i),'color',[0.5 0.5 0.5],'HorizontalAlignment', 'center');
        h.Bearing{i}.Line  = plot3(P.Bearing{i}.z*[1 1], -[rPlot+0.01 rText-0.025], [0 0],'color',[0.5 0.5 0.5]);
    end

    %and the nodes
    for i = length(P.Rotor):-1:1
        rText = rText + 0.05;
        for j = 1:length(P.Rotor{i}.Nodes)
            h.Rotor{i}.Nodes(j).Label = text(P.Rotor{i}.Nodes(j),rText,0,sprintf('Node %d',j),'color',col(i,:),'HorizontalAlignment', 'center');
            h.Rotor{i}.Nodes(j).Line = plot3(P.Rotor{i}.Nodes(j)*[1 1], [rNodes{i}(j)+0.01 rText-0.025], [0 0],'color',col(i,:));
        end
        hLeg(i) = h.Rotor{i}.Nodes(1).Line;
    end

    h.Legend = legend(hLeg,Rotor_Names,'Location','EastOutside');
end

xlabel('z (m)')
ylabel('x (m)')
zlabel('y (m)')
view(2)
axis equal

if nargout > 0
    varargout{1} = h;
end

function [X,Y,Z] = plot_cylinder(ri,ro,z0,t)
r = [ri ro ro ri];
z = [-1 -1 1 1]*t/2 + z0;

%go back to start if hollow
if ri > 0
    r(end+1) = r(1);
    z(end+1) = z(1);
end

[X,Y,~] = cylinder(r,16);
Z = repmat(z(:),1,size(X,2));