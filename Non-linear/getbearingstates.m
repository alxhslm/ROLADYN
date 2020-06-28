function bearing_states = getbearingstates(States,P,hbm)
NPts = size(States.x,2);

bearing_states.x     = P.Model.Bearing.S*States.x(1:P.Model.NDof,:);
bearing_states.xdot  = P.Model.Bearing.S*States.xdot(1:P.Model.NDof,:);
bearing_states.xddot = P.Model.Bearing.S*States.xddot(1:P.Model.NDof,:);

bearing_states.u     = States.u;
bearing_states.udot  = States.udot;
bearing_states.uddot = States.uddot;

xB     = States.x(P.Model.NDof+1:end,:);
xdotB  = States.xdot(P.Model.NDof+1:end,:);
xddotB = States.xddot(P.Model.NDof+1:end,:);

bearing_states.xInt     = zeros(P.Model.NDofInt,NPts);
bearing_states.xdotInt  = zeros(P.Model.NDofInt,NPts);
bearing_states.xddotInt = zeros(P.Model.NDofInt,NPts);

iDofIn = 0;
iDofOut = 0;
for i = 1:length(P.Bearing)
    if strcmp(P.Bearing{i}.Model,'REB') && P.Bearing{i}.Params.Model.NDof > 0
        REB = P.Bearing{i}.Params;

        if P.Model.bCompressREB && (size(xB,1) == P.Model.Reduced.NDofInt)
            [xInt,xdotInt,xddotInt] = shift_balls_fast(xB(iDofIn+(1:REB.Model.NDof),:),xdotB(iDofIn+(1:REB.Model.NDof),:),xddotB(iDofIn+(1:REB.Model.NDof),:),REB,P.Model.iRot,hbm);
            iDofIn  = iDofIn  + REB.Model.NDof;
        else
            xInt = xB(iDofIn+(1:P.Bearing{i}.NDofInt),:);
            xdotInt(iDofOut+(1:P.Bearing{i}.NDofInt),:)  =  xdotB(iDofIn+(1:P.Bearing{i}.NDofInt),:);
            xddotInt(iDofOut+(1:P.Bearing{i}.NDofInt),:) =  xddotB(iDofIn+(1:P.Bearing{i}.NDofInt),:);
            iDofIn  = iDofIn  + P.Bearing{i}.NDofInt;
        end

        bearing_states.xInt(iDofOut+(1:P.Bearing{i}.NDofInt),:)     = xInt;
        bearing_states.xdotInt(iDofOut+(1:P.Bearing{i}.NDofInt),:)  = xdotInt;
        bearing_states.xddotInt(iDofOut+(1:P.Bearing{i}.NDofInt),:) = xddotInt;

        iDofOut = iDofOut + P.Bearing{i}.NDofInt;

    else

        bearing_states.xInt(iDofOut+(1:P.Bearing{i}.NDofInt),:)     =  xB(iDofIn+(1:P.Bearing{i}.NDofInt),:);
        bearing_states.xdotInt(iDofOut+(1:P.Bearing{i}.NDofInt),:)  =  xdotB(iDofIn+(1:P.Bearing{i}.NDofInt),:);
        bearing_states.xddotInt(iDofOut+(1:P.Bearing{i}.NDofInt),:) =  xddotB(iDofIn+(1:P.Bearing{i}.NDofInt),:);

        iDofIn  = iDofIn  + P.Bearing{i}.NDofInt;
        iDofOut = iDofOut + P.Bearing{i}.NDofInt;
    end
end

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
