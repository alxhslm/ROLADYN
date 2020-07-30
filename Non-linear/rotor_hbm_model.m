function varargout = rotor_hbm_model(part,States,hbm,problem)
O  = States.w0(1);
P = problem.P;

bearing_states = getbearingstates(States,P,hbm);
bearing_states.A = O*States.t(P.Model.iRot,:);
bearing_states.O = O + 0*States.t(P.Model.iRot,:);
bearing_states.bSolve = 0;
bearing_states.bNL = 1;
bearing_states.bLin = 0;

States.xInt = States.x(P.Model.NDof + 1:end,:);
States.x = States.x(1:P.Model.NDof,:);
States.xdot = States.xdot(1:P.Model.NDof,:);
States.xddot = States.xddot(1:P.Model.NDof,:);

varargout = {};
switch part
    case  'nl'        
        if P.Model.bNL
            Forces = bearingforces(P,bearing_states);
        else
            Forces = linear_bearingforces(P,bearing_states);
        end
        
        Fi = Forces.FInt(P.Model.Reduced.iInt,:);
        Fnl = P.Model.Bearing.S'*Forces.F;
        
        varargout{end+1} = [Fnl; Fi];
    case 'output'
        if P.Model.bNL
            Forces = bearingforces(P,bearing_states);
        else
            Forces = linear_bearingforces(P,bearing_states);
        end
        varargout{end+1} = Forces.F;
    otherwise
        if P.Model.bAnalyticalDerivs
            % ///// Analytical derivs /////
            if P.Model.bNL
                [~,Stiffness] = bearingforces(P,bearing_states);
            else
                [~,Stiffness] = linear_bearingforces(P,bearing_states);
            end
        else
            f0 = rotor_hbm_model('nl',States,hbm,problem);
        end
        switch part
            case 'nl_x'
                if P.Model.bAnalyticalDerivs
                    Kqq =  mtimesx(P.Model.Bearing.S',mtimesx(Stiffness.Kqq,P.Model.Bearing.S));
                    Kqx =  mtimesx(P.Model.Bearing.S',Stiffness.Kqx);
                    Kxq =  mtimesx(Stiffness.Kxq,P.Model.Bearing.S);
                    Kxx =  Stiffness.Kxx;
                    K = [Kqq Kqx; Kxq Kxx];
                else
                    K = numerical_deriv(States,f0,'x',hbm,problem);
                end
                varargout{end+1} = K;
            case 'nl_xdot'
                if P.Model.bAnalyticalDerivs
                    Cqq =  mtimesx(P.Model.Bearing.S',mtimesx(Stiffness.Cqq,P.Model.Bearing.S));
                    Cqx =  mtimesx(P.Model.Bearing.S',Stiffness.Cqx);
                    Cxq =  mtimesx(Stiffness.Cxq,P.Model.Bearing.S);
                    Cxx =  Stiffness.Cxx;
                    C = [Cqq Cqx; Cxq Cxx];
                else
                    C = numerical_deriv(States,f0,'xdot',hbm,problem);
                end
                varargout{end+1} = C;
            case 'nl_xddot'
                if P.Model.bAnalyticalDerivs
                    Mqq =  mtimesx(P.Model.Bearing.S',mtimesx(Stiffness.Mqq,P.Model.Bearing.S));
                    Mqx =  mtimesx(P.Model.Bearing.S',Stiffness.Mqx);
                    Mxq =  mtimesx(Stiffness.Mxq,P.Model.Bearing.S);
                    Mxx =  Stiffness.Mxx;
                    M = [Mqq Mqx; Mxq Mxx];
                else
                    M = numerical_deriv(States,f0,'xddot',hbm,problem);
                end
                varargout{end+1} = M;
            case 'nl_u'
                if P.Model.bAnalyticalDerivs
                    Kqu =  mtimesx(P.Model.Bearing.S',Stiffness.Kqu);
                    Kxu =  Stiffness.Kxu;
                    Ku = [Kqu; Kxu];
                else
                    Ku = numerical_deriv(States,f0,'u',hbm,problem);
                end
                varargout{end+1} = Ku;
            case 'nl_udot'
                if P.Model.bAnalyticalDerivs
                    Cqu =  mtimesx(P.Model.Bearing.S',Stiffness.Cqu);
                    Cxu =  Stiffness.Cxu;
                    Cu = [Cqu; Cxu];
                else
                    Cu = numerical_deriv(States,f0,'udot',hbm,problem);
                end
                varargout{end+1} = Cu;
            case 'nl_uddot'
                if P.Model.bAnalyticalDerivs
                    Mqu =  mtimesx(P.Model.Bearing.S',Stiffness.Mqu);
                    Mxu =  Stiffness.Mxu;
                    Mu = [Mqu; Mxu];
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
           P.Mesh.Bearing.Kqq*(bearing_states.x-P.Mesh.Bearing.xb) + P.Mesh.Bearing.Cqq*bearing_states.xdot + P.Mesh.Bearing.Mqq*bearing_states.xddot + ...
           P.Mesh.Bearing.Kqx*(bearing_states.xInt-P.Mesh.xInt) + P.Mesh.Bearing.Cqx*bearing_states.xdotInt + P.Mesh.Bearing.Mqx*bearing_states.xddotInt + ...
           P.Mesh.Bearing.Kqu*bearing_states.u + P.Mesh.Bearing.Cqu*bearing_states.udot + P.Mesh.Bearing.Mqu*bearing_states.uddot;

Forces.FInt = P.Mesh.Bearing.Kxq*(bearing_states.x-P.Mesh.Bearing.xb) + P.Mesh.Bearing.Cxq*bearing_states.xdot + P.Mesh.Bearing.Mxq*bearing_states.xddot + ...
           P.Mesh.Bearing.Kxx*(bearing_states.xInt-P.Mesh.xInt) + P.Mesh.Bearing.Cxx*bearing_states.xdotInt + P.Mesh.Bearing.Mxx*bearing_states.xddotInt + ...
           P.Mesh.Bearing.Kxu*bearing_states.u + P.Mesh.Bearing.Cxu*bearing_states.udot + P.Mesh.Bearing.Mxu*bearing_states.uddot;

            
if nargout > 1
    Stiffness.Kqq = repmat(P.Mesh.Bearing.Kqq,1,1,NPts);
    Stiffness.Kqx = repmat(P.Mesh.Bearing.Kqx,1,1,NPts);
    Stiffness.Kxq = repmat(P.Mesh.Bearing.Kxq,1,1,NPts);
    Stiffness.Kxx = repmat(P.Mesh.Bearing.Kxx,1,1,NPts);

    Stiffness.Cqq = repmat(P.Mesh.Bearing.Cqq,1,1,NPts);
    Stiffness.Cqx = repmat(P.Mesh.Bearing.Cqx,1,1,NPts);
    Stiffness.Cxq = repmat(P.Mesh.Bearing.Cxq,1,1,NPts);
    Stiffness.Cxx = repmat(P.Mesh.Bearing.Cxx,1,1,NPts);

    Stiffness.Mqq = repmat(P.Mesh.Bearing.Mqq,1,1,NPts);
    Stiffness.Mqx = repmat(P.Mesh.Bearing.Mqx,1,1,NPts);
    Stiffness.Mxq = repmat(P.Mesh.Bearing.Mxq,1,1,NPts);
    Stiffness.Mxx = repmat(P.Mesh.Bearing.Mxx,1,1,NPts);

    Stiffness.Kqu = repmat(P.Mesh.Bearing.Kqu,1,1,NPts);
    Stiffness.Kxu = repmat(P.Mesh.Bearing.Kxu,1,1,NPts);

    Stiffness.Cqu = repmat(P.Mesh.Bearing.Cqu,1,1,NPts);
    Stiffness.Cxu = repmat(P.Mesh.Bearing.Cxu,1,1,NPts);

    Stiffness.Mqu = repmat(P.Mesh.Bearing.Mqu,1,1,NPts);
    Stiffness.Mxu = repmat(P.Mesh.Bearing.Mxu,1,1,NPts);
end