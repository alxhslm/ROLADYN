function [Fc,Fi,Fo,Mg] = dynamic_ball_loads(B,alpha_i,alpha_o,Oi,Oo)

if B.Options.bComplexDynLoads
    %select the ball speed ratio depending on the operating condition
    % alpha_i = angle_fix(alpha_i);
    % alpha_o = angle_fix(alpha_o);
    % alpha_d = wrapTo2Pi(alpha_i - alpha_o);
    alpha_d = alpha_i - alpha_o;
    switch B.Options.Control
        case 'outer'
            rCagei = (1 - B.Geometry.gamma * cos(alpha_i))./(1+cos(alpha_d));
            beta  = atan(sin(alpha_o)./(cos(alpha_o)+B.Geometry.gamma));
            rRoll = 1./((cos(alpha_o)+tan(beta).*sin(alpha_o))./((1+B.Geometry.gamma*cos(alpha_o) + (cos(alpha_i)+tan(beta).*sin(alpha_i))./(1-B.Geometry.gamma*cos(alpha_i))).*(B.Geometry.gamma*cos(beta))));
        case 'inner'
            rCagei = (cos(alpha_d) - B.Geometry.gamma * cos(alpha_i))./(1+cos(alpha_d));
            beta = atan(sin(alpha_i)./(cos(alpha_i)-B.Geometry.gamma));
            rRoll = -1./((cos(alpha_o)+tan(beta).*sin(alpha_o))./((1+B.Geometry.gamma*cos(alpha_o) + (cos(alpha_i)+tan(beta).*sin(alpha_i))./(1-B.Geometry.gamma*cos(alpha_i))).*(B.Geometry.gamma*cos(beta))));
        otherwise
            error('Unknown raceway condition')
    end
else
    sgn =  sign(alpha_i+eps);
    wons =  0*alpha_i+1;
    rCagei = B.Kinematics.rCagei*wons;
    switch B.Options.Control
        case 'outer'
            beta = B.Kinematics.betao*sgn;
            rRoll = B.Kinematics.rRollo*wons;
        case 'inner'
            beta = B.Kinematics.betai*sgn;
            rRoll = B.Kinematics.rRolli*wons;
    end
end

%compute the roll and cage speeds
Ocage = rCagei.*Oi + (1-rCagei).*Oo;
wRoll = rRoll .* (Oo - Oi);
wdotPivot = 0*alpha_i;

%compute centrifugal, gyroscopic and pivotal loads
Fc = B.Options.bCentrifugal*B.Dynamics.mb*(B.Geometry.dm/2)*Ocage.^2;
Mg = B.Options.bGyro*B.Dynamics.Jb*wRoll.*Ocage.*sin(beta);
Mp = B.Options.bPivot*B.Dynamics.Jb*wdotPivot;

%split moment reaction forces between inner/outer contacts
if ~B.Options.bGyro
    Fi = 0*alpha_i;
    Fo = 0*alpha_i;
else
    switch B.Options.Control
        case 'outer'
            Fo = (Mg-Mp)*2/B.Geometry.D;
            Fi = 0*alpha_i;
        case 'inner'
            Fo = 0*alpha_i;
            Fi = (Mg-Mp)*2/B.Geometry.D;
    end
end