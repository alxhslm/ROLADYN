function [usi,uso] = ehd_sum_speeds(B,alpha_i,alpha_o,Oi,Oo)

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
            rRoll = 1/((cos(alpha_o)+tan(beta).*sin(alpha_o))./((1+B.Geometry.gamma*cos(alpha_o) + (cos(alpha_i)+tan(beta).*sin(alpha_i))./(1-B.Geometry.gamma*cos(alpha_i))).*(B.Geometry.gamma*cos(beta))));
            %             rPivot = 1;
        case 'inner'
            rCagei = (cos(alpha_d) - B.Geometry.gamma * cos(alpha_i))./(1+cos(alpha_d));
            beta = atan(sin(alpha_i)./(cos(alpha_i)-B.Geometry.gamma));
            rRoll = -1/((cos(alpha_o)+tan(beta).*sin(alpha_o))./((1+B.Geometry.gamma*cos(alpha_o) + (cos(alpha_i)+tan(beta).*sin(alpha_i))./(1-B.Geometry.gamma*cos(alpha_i))).*(B.Geometry.gamma*cos(beta))));
            %             rPivot = 1;
        otherwise
            error('Unknown raceway condition')
    end
else
    rCagei = B.Kinematics.rCagei;    
    switch B.Options.Control
        case 'outer'
            rRoll = B.Kinematics.rRollo;
        case 'inner'
            rRoll = B.Kinematics.rRolli;
    end
end
    %compute the roll and cage speeds
    Ocage = rCagei.*Oi + (1-rCagei).*Oo;
    wRoll = rRoll .* (Oo - Oi);    

usi = (Oi - Ocage)*B.Geometry.di - (wRoll - Ocage)*B.Geometry.D/2;
uso = (Oi - Ocage)*B.Geometry.do + (wRoll - Ocage)*B.Geometry.D/2;