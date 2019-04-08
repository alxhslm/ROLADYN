function varargout = rotor_hbm_model(part,States,hbm,problem)
O  = States.w0(1)/hbm.harm.rFreqRatio(1);
P = problem.P;

NPts = size(States.x,2);

bearing_states = getbearingstates(States,O,P,hbm);
bearing_states.A = O*States.t(P.Model.iRot,:);
bearing_states.O = O + 0*States.t(P.Model.iRot,:);
bearing_states.bSolve = 0;

varargout = {};
switch part
    case  {'nl','all','output'}
        Fgyro = O*P.Model.Rotor.G*States.x(1:P.Model.NDof,:);
        Forces = bearingforces(P,bearing_states);
        
        if P.Model.bNL
            Fb = P.Model.A'*Forces.F - P.Model.Bearing.C*States.xdot(1:P.Model.NDof,:);
        else
            Fb = P.Model.Bearing.F0 + P.Model.Bearing.K*(States.x(1:P.Model.NDof,:)-P.Model.x0);% + P.Model.Bearing.C*xdotCG;
        end
        
%         Fnl = P.Model.A'*Forces.F;
%         Flin = P.Model.Bearing.F0 + P.Model.Bearing.K*(xCG-P.Model.x0) + P.Model.Bearing.C*xdotCG;
        
        Fnl = Fb + Fgyro;
        if P.Model.bCompressREB
            Fi = Forces.FInt(P.Model.Reduced.iInt,:);
        else
            Fi = Forces.FInt;
        end
        
        if strcmp(part,'output')
            varargout{end+1} = P.Model.A'*Forces.F;
        elseif strcmp(part,'nl')
                varargout{end+1} = [Fnl; Fi];
       
        elseif strcmp(part,'all')
             varargout{end+1} = [Fnl; Fi];
        else
            error('Unknown option');
        end
    otherwise
        Cgyro = O*P.Model.Rotor.G;
        if P.Model.bAnalyticalDerivs
            [Forces,Stiffness] = bearingforces(P,bearing_states);
            
            Kqq =  mtimesx(P.Model.A',mtimesx(Stiffness.Kqq,P.Model.A));
            Kqx =  mtimesx(P.Model.A',Stiffness.Kqx);
            Kxq =  mtimesx(Stiffness.Kxq,P.Model.A);
            Kxx =  Stiffness.Kxx;
                
            Cqq =  mtimesx(P.Model.A',mtimesx(Stiffness.Cqq,P.Model.A)) + Cgyro;
            Cqx =  mtimesx(P.Model.A',Stiffness.Cqx);
            Cxq =  mtimesx(Stiffness.Cxq,P.Model.A);
            Cxx =  Stiffness.Cxx;
            
            Kqu = mtimesx(P.Model.A',mtimesx(Stiffness.Kqu,P.Mesh.Excite.Sgd));
            Kxu = mtimesx(Stiffness.Kxu,P.Mesh.Excite.Sgd);
            
            Cqu = mtimesx(P.Model.A',mtimesx(Stiffness.Cqu,P.Mesh.Excite.Sgd));
            Cxu = mtimesx(Stiffness.Cxu,P.Mesh.Excite.Sgd);
        else
            f0 = rotor_hbm_model('all',States,hbm,problem);
            h = 1E-10;
            
            Kqq = zeros(P.Model.NDof,P.Model.NDof,NPts);
            Kxq = zeros(P.Model.NDofInt,P.Model.NDof,NPts);
            Cqq = zeros(P.Model.NDof,P.Model.NDof,NPts);
            Cxq = zeros(P.Model.NDofInt,P.Model.NDof,NPts);
            x0 = States.x;
            xdot0 = States.xdot;
            for i = 1:P.Model.NDof
                States.x(i,:) = States.x(i,:) + h;
                f = rotor_hbm_model('all',States,hbm,problem);
                Kqq(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Kxq(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                States.x = x0;
                
                States.xdot(i,:) = States.xdot(i,:) + h;
                f = rotor_hbm_model('all',States,hbm,problem);
                Cqq(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Cxq(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                States.xdot = xdot0;
            end

            Kqx = zeros(P.Model.NDof,P.Model.NDofInt,NPts);
            Kxx = zeros(P.Model.NDofInt,P.Model.NDofInt,NPts);
            Cqx = zeros(P.Model.NDof,P.Model.NDofInt,NPts);
            Cxx = zeros(P.Model.NDofInt,P.Model.NDofInt,NPts);
            
            Kqu = zeros(P.Model.NDof,P.Mesh.NExcite,NPts);
            Kxu = zeros(P.Model.NDofInt,P.Mesh.NExcite,NPts);
            Cqu = zeros(P.Model.NDof,P.Mesh.NExcite,NPts);
            Cxu = zeros(P.Model.NDofInt,P.Mesh.NExcite,NPts);
            u0 = States.u;
            udot0 = States.udot;
            for i = 1:P.Mesh.NExcite
                States.u(i,:) = States.u(i,:) + h;
                f = rotor_hbm_model('all',States,hbm,problem);
                Kqu(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Kxu(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                States.u = u0;
                
                States.udot(i,:) = States.udot(i,:) + h;
                f = rotor_hbm_model('all',States,hbm,problem);
                Cqu(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Cxu(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                States.udot = udot0;
            end
        end
        
        switch part
            case 'nl_x'
                varargout{end+1} = [Kqq Kqx; Kxq Kxx];
            case 'nl_xdot'
                varargout{end+1} = [Cqq Cqx; Cxq Cxx];
            case 'nl_u'
                varargout{end+1} = [Kqu; Kxu];         
            case 'nl_udot'
                varargout{end+1} = [Cqu; Cxu]; 
            case 'nl_xddot'
                varargout{end+1} = 0*[Kqx; Kxx];
            case 'nl_uddot'
                varargout{end+1} = 0*[Kqu; Kxu];
            case 'nl_w'
                varargout{end+1} = [];
                %too complicated to do analytically
        end
end