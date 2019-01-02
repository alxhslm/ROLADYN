function States = getbearingstates(x,xdot,xddot,u,udot,uddot,xalg,O,P,hbm)
States.x     = P.Model.A*x;
States.xdot  = P.Model.A*xdot;
States.xddot = P.Model.A*xddot;

xB     = xalg;
xdotB  = 0*xalg;
xddotB = 0*xalg;

States.xInt     = zeros(P.Model.NDofInt,size(x,2));
States.xdotInt  = zeros(P.Model.NDofInt,size(x,2));
States.xddotInt = zeros(P.Model.NDofInt,size(x,2));

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
                                    
            States.xInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:)     = xInt;
            States.xdotInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:)  = xdotInt;
            States.xddotInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:) = xddotInt;
            
            iDofOut = iDofOut + P.Bearing{i}.NDofInt(j);
            
        else
            
            States.xInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:)     =  xB(iDofIn+(1:P.Bearing{i}.NDofInt(j)),:);
            States.xdotInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:)  =  xdotB(iDofIn+(1:P.Bearing{i}.NDofInt(j)),:);
            States.xddotInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:) =  xddotB(iDofIn+(1:P.Bearing{i}.NDofInt(j)),:);
            
            iDofIn  = iDofIn  + P.Bearing{i}.NDofInt(j);
            iDofOut = iDofOut + P.Bearing{i}.NDofInt(j);
        end
    end
end

States.u     = P.Mesh.Excite.Sgd*u;
States.udot  = P.Mesh.Excite.Sgd*udot;
States.uddot = P.Mesh.Excite.Sgd*uddot;


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