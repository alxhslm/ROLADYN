function bearing_states = getbearingstates(States,O,P,hbm)
NPts = size(States.x,2);

bearing_states.x     = P.Model.A*States.x(1:P.Model.NDof,:);
bearing_states.xdot  = P.Model.A*States.xdot(1:P.Model.NDof,:);
bearing_states.xddot = P.Model.A*States.xddot(1:P.Model.NDof,:);

xB     = States.x(P.Model.NDof+1:end,:);
xdotB  = 0*xB;
xddotB = 0*xB;

bearing_states.xInt     = zeros(P.Model.NDofInt,NPts);
bearing_states.xdotInt  = zeros(P.Model.NDofInt,NPts);
bearing_states.xddotInt = zeros(P.Model.NDofInt,NPts);

iDofIn = 0;
iDofOut = 0;
for i = 1:length(P.Bearing)
    for j = 1:2
        if strcmp(P.Bearing{i}.Model{j},'REB') && P.Bearing{i}.Params{j}.Model.NDof > 0
            REB = P.Bearing{i}.Params{j};
           
            if P.Model.bCompressREB && (size(xB,1) == P.Model.Reduced.NDofInt)
%                 [xInt,xdotInt,xddotInt] = shift_balls(xB(iDofIn+(1:REB.Model.NDof),:),xdotB(iDofIn+(1:REB.Model.NDof),:),xddotB(iDofIn+(1:REB.Model.NDof),:),REB,P.Model.iRot,hbm);
                [xInt,xdotInt,xddotInt] = shift_balls_fast(xB(iDofIn+(1:REB.Model.NDof),:),xdotB(iDofIn+(1:REB.Model.NDof),:),xddotB(iDofIn+(1:REB.Model.NDof),:),REB,P.Model.iRot,hbm);
                iDofIn  = iDofIn  + REB.Model.NDof;
            else
                xInt = xB(iDofIn+(1:P.Bearing{i}.NDofInt(j)),:);
                xdotInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:)  =  xdotB(iDofIn+(1:P.Bearing{i}.NDofInt(j)),:);
                xddotInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:) =  xddotB(iDofIn+(1:P.Bearing{i}.NDofInt(j)),:);
                iDofIn  = iDofIn  + P.Bearing{i}.NDofInt(j);
            end
                                    
            bearing_states.xInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:)     = xInt;
            bearing_states.xdotInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:)  = xdotInt;
            bearing_states.xddotInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:) = xddotInt;
            
            iDofOut = iDofOut + P.Bearing{i}.NDofInt(j);
            
        else
            
            bearing_states.xInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:)     =  xB(iDofIn+(1:P.Bearing{i}.NDofInt(j)),:);
            bearing_states.xdotInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:)  =  xdotB(iDofIn+(1:P.Bearing{i}.NDofInt(j)),:);
            bearing_states.xddotInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:) =  xddotB(iDofIn+(1:P.Bearing{i}.NDofInt(j)),:);
            
            iDofIn  = iDofIn  + P.Bearing{i}.NDofInt(j);
            iDofOut = iDofOut + P.Bearing{i}.NDofInt(j);
        end
    end
end

bearing_states.u     = P.Mesh.Excite.Sgd*States.u;
bearing_states.udot  = P.Mesh.Excite.Sgd*States.udot;
bearing_states.uddot = P.Mesh.Excite.Sgd*States.uddot;


function [xInt,xdotInt,xddotInt] = shift_balls_fast(x,xdot,xddot,REB,iRot,hbm)
xInt     = zeros(REB.Model.NDofTot,prod(hbm.harm.Nfft));
xdotInt  = zeros(REB.Model.NDofTot,prod(hbm.harm.Nfft));
xddotInt = zeros(REB.Model.NDofTot,prod(hbm.harm.Nfft));

Z = REB.Setup.Z;
for i = 1:Z
    if iRot == 1
       kShift = -(i-1)/Z * hbm.harm.Nfft(1); 
    elseif iRot == 2
        kShift = -(i-1)/Z * hbm.harm.Nfft(2) *  hbm.harm.Nfft(1);
    end

    xInt(i + ((1:REB.Model.NDof)-1)*Z,:) = circshift(x,    kShift,2);
    xdotInt(i + ((1:REB.Model.NDof)-1)*Z ,:) = circshift(xdot, kShift,2);
    xddotInt(i + ((1:REB.Model.NDof)-1)*Z,:) = circshift(xddot,kShift,2);
end

function [xInt,xdotInt,xddotInt] = shift_balls(x,xdot,xddot,REB,iRot,hbm)
x     = reshape(x'    ,hbm.harm.Nfft(1),hbm.harm.Nfft(2),REB.Model.NDof);
xdot  = reshape(xdot' ,hbm.harm.Nfft(1),hbm.harm.Nfft(2),REB.Model.NDof);
xddot = reshape(xddot',hbm.harm.Nfft(1),hbm.harm.Nfft(2),REB.Model.NDof);

xB     = zeros(hbm.harm.Nfft(1),hbm.harm.Nfft(2),REB.Model.NDofTot);
xdotB  = zeros(hbm.harm.Nfft(1),hbm.harm.Nfft(2),REB.Model.NDofTot);
xddotB = zeros(hbm.harm.Nfft(1),hbm.harm.Nfft(2),REB.Model.NDofTot);

Z = REB.Setup.Z;

for i = 1:Z
    kShift = -(i-1)/Z * hbm.harm.Nfft(iRot);
    for j = 1:REB.Model.NDof
        xB(:,:,    i + (j-1)*Z) = circshift(x(:,:,j),    kShift,iRot);
        xdotB(:,:, i + (j-1)*Z) = circshift(xdot(:,:,j), kShift,iRot);
        xddotB(:,:,i + (j-1)*Z) = circshift(xddot(:,:,j),kShift,iRot);
    end
end
xInt     = reshape(xB,    prod(hbm.harm.Nfft),REB.Model.NDofTot)';
xdotInt  = reshape(xdotB, prod(hbm.harm.Nfft),REB.Model.NDofTot)';
xddotInt = reshape(xddotB,prod(hbm.harm.Nfft),REB.Model.NDofTot)';