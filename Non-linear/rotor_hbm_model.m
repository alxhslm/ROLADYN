function varargout = rotor_hbm_model(part,States,hbm,problem)
O  = States.w0(1);
P = problem.P;

bearing_states.x = P.Model.Bearing.S*States.x;
bearing_states.xdot = P.Model.Bearing.S*States.xdot;
bearing_states.xddot = P.Model.Bearing.S*States.xddot;

bearing_states.A = O*States.t(P.Model.iRot,:);
bearing_states.O = O + 0*States.t(P.Model.iRot,:);

varargout = {};
switch part
    case  'nl'        
        if P.Model.bNL
            Forces = bearingforces(P,bearing_states, 1, 0);
        else
            Forces = linear_bearingforces(P,bearing_states);
        end
        
        Fnl = P.Model.Bearing.S'*Forces.F;
        
        varargout{end+1} = Fnl;
    case 'output'
        if P.Model.bNL
            Forces = bearingforces(P,bearing_states, 1, 0);
        else
            Forces = linear_bearingforces(P,bearing_states);
        end
        varargout{end+1} = Forces.F;
    otherwise
        if P.Model.bAnalyticalDerivs
            % ///// Analytical derivs /////
            if P.Model.bNL
                [~,Stiffness] = bearingforces(P,bearing_states, 1, 0);
            else
                [~,Stiffness] = linear_bearingforces(P,bearing_states);
            end
        else
            f0 = rotor_hbm_model('nl',States,hbm,problem);
        end
        switch part
            case 'nl_x'
                if P.Model.bAnalyticalDerivs
                    K =  mtimesx(P.Model.Bearing.S',mtimesx(Stiffness.K,P.Model.Bearing.S));
                else
                    K = numerical_deriv(States,f0,'x',hbm,problem);
                end
                varargout{end+1} = K;
            case 'nl_xdot'
                if P.Model.bAnalyticalDerivs
                    C =  mtimesx(P.Model.Bearing.S',mtimesx(Stiffness.C,P.Model.Bearing.S));
                else
                    C = numerical_deriv(States,f0,'xdot',hbm,problem);
                end
                varargout{end+1} = C;
            case 'nl_xddot'
                if P.Model.bAnalyticalDerivs
                    M =  mtimesx(P.Model.Bearing.S',mtimesx(Stiffness.M,P.Model.Bearing.S));
                else
                    M = numerical_deriv(States,f0,'xddot',hbm,problem);
                end
                varargout{end+1} = M;
            case 'nl_u'
                if P.Model.bAnalyticalDerivs
                    Ku =  mtimesx(P.Model.Bearing.S',Stiffness.Ku);
                else
                    Ku = numerical_deriv(States,f0,'u',hbm,problem);
                end
                varargout{end+1} = Ku;
            case 'nl_udot'
                if P.Model.bAnalyticalDerivs
                    Cu =  mtimesx(P.Model.Bearing.S',Stiffness.Cu);
                else
                    Cu = numerical_deriv(States,f0,'udot',hbm,problem);
                end
                varargout{end+1} = Cu;
            case 'nl_uddot'
                if P.Model.bAnalyticalDerivs
                    Mu =  mtimesx(P.Model.Bearing.S',Stiffness.Mu);
                else
                    Mu = numerical_deriv(States,f0,'uddot',hbm,problem);
                end
                varargout{end+1} = Mu;
            case 'nl_w'
                if P.Model.bAnalyticalDerivs
                    %only valid if VC is turned off
                    Dw = repmat({0*States.x},1,2);
                else
                    w0 = States.w0;
                    h = 1E-6;
                    Dw = cell(1,2);
                    for i = 1:length(w0)
                        States.w0(i) = States.w0(i) + h;
                        Dw{i} = (rotor_hbm_model('nl',States,hbm,problem) - f0)/h;
                        States.w0 = w0;
                    end
                end
                varargout{end+1} = Dw;
        end
end

function K = numerical_deriv(States,f0,field,hbm,problem)
NPts = size(States.(field),2);
NInput = size(States.(field),1);
K = zeros(problem.NDof,NInput,NPts);
if hbm.dependence.(field)
    x0 = States.(field);
    h = 1E-10;
    for i = 1:NInput
        States.(field)(i,:) = States.(field)(i,:) + h;
        f = rotor_hbm_model('nl',States,hbm,problem);
        K(:,i,:) = permute((f-f0)/h,[1 3 2]);
        States.(field) = x0;
    end
end

function [Forces,Stiffness] = linear_bearingforces(P,bearing_states)
NPts = size(bearing_states.x,2);

Forces.F = P.Mesh.Bearing.Fb + ...
           P.Mesh.Bearing.K*(bearing_states.x-P.Mesh.Bearing.xb) + P.Mesh.Bearing.C*bearing_states.xdot + P.Mesh.Bearing.M*bearing_states.xddot + ...
           P.Mesh.Bearing.Cu*bearing_states.u + P.Mesh.Bearing.Cu*bearing_states.udot + P.Mesh.Bearing.Mu*bearing_states.uddot;
            
if nargout > 1
    Stiffness.K = repmat(P.Mesh.Bearing.K,1,1,NPts);
    Stiffness.C = repmat(P.Mesh.Bearing.C,1,1,NPts);
    Stiffness.M = repmat(P.Mesh.Bearing.M,1,1,NPts);

    Stiffness.Ku = repmat(P.Mesh.Bearing.Ku,1,1,NPts);
    Stiffness.Cu = repmat(P.Mesh.Bearing.Cu,1,1,NPts);
    Stiffness.Mu = repmat(P.Mesh.Bearing.Mu,1,1,NPts);
end