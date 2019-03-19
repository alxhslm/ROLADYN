function varargout = rotor_hbm_model(part,t,x,xdot,xddot,u,udot,uddot,hbm,problem,w0)
P = problem.P;
O  = w0(1)/hbm.harm.rFreqRatio(1);
NPts = size(x,2);
tRot = t(P.Model.iRot,:);

NDof = size(x,1);
NInput = size(u,1);


xInt = x(P.Model.NDof+1:end,:);
xCG = x(1:P.Model.NDof,:);
xdotCG = xdot(1:P.Model.NDof,:);
xddotCG = xddot(1:P.Model.NDof,:);

if size(xInt,1) == P.Model.NDofInt
    bReduced = 0;
elseif size(xInt,1) == P.Model.Reduced.NDofInt && P.Model.bCompressREB
    bReduced = 1;
else
    error('Wrong size')
end

States = getbearingstates(xCG,xdotCG,xddotCG,u,udot,uddot,xInt,O,P,hbm);
States.A = O*tRot;
States.O = O + 0*tRot;
States.bSolve = 0;

x = [xCG; States.xInt];
xdot = [xdotCG; States.xdotInt];
xddot = [xddotCG; States.xddotInt];

varargout = {};
switch part
    case  {'nl','all','output'}
        Fgyro = O*P.Model.Rotor.G*xCG;
        Forces = bearingforces(P,States);
        
        if P.Model.bNL
            Fb = P.Model.A'*Forces.F - P.Model.Bearing.C*xdotCG;
        else
            Fb = P.Model.Bearing.F0 + P.Model.Bearing.K*(xCG-P.Model.x0);% + P.Model.Bearing.C*xdotCG;
        end
        
%         Fnl = P.Model.A'*Forces.F;
%         Flin = P.Model.Bearing.F0 + P.Model.Bearing.K*(xCG-P.Model.x0) + P.Model.Bearing.C*xdotCG;
        
        Fnl = Fb + Fgyro;
        if bReduced
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
        if isfield(P.Model,'bAnalyticalDerivs') && P.Model.bAnalyticalDerivs
            [Forces,Stiffness] = bearingforces(P,States);
            
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
            f0 = rotor_hbm_model('all',t,x,xdot,xddot,u,udot,uddot,hbm,problem,w0);
            h = 1E-10;
            
            Kqq = zeros(P.Model.NDof,P.Model.NDof,NPts);
            Kxq = zeros(P.Model.NDofInt,P.Model.NDof,NPts);
            Cqq = zeros(P.Model.NDof,P.Model.NDof,NPts);
            Cxq = zeros(P.Model.NDofInt,P.Model.NDof,NPts);
            x0 = x;
            xdot0 = xdot;
            for i = 1:P.Model.NDof
                x(i,:) = x(i,:) + h;
                f = rotor_hbm_model('all',t,x,xdot,xddot,u,udot,uddot,hbm,problem,w0);
                Kqq(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Kxq(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                x = x0;
                
                xdot(i,:) = xdot(i,:) + h;
                f = rotor_hbm_model('all',t,x,xdot,xddot,u,udot,uddot,hbm,problem,w0);
                Cqq(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Cxq(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                xdot = xdot0;
            end

            Kqx = zeros(P.Model.NDof,P.Model.NDofInt,NPts);
            Kxx = zeros(P.Model.NDofInt,P.Model.NDofInt,NPts);
            Cqx = zeros(P.Model.NDof,P.Model.NDofInt,NPts);
            Cxx = zeros(P.Model.NDofInt,P.Model.NDofInt,NPts);
            
            Kqu = zeros(P.Model.NDof,P.Mesh.NExcite,NPts);
            Kxu = zeros(P.Model.NDofInt,P.Mesh.NExcite,NPts);
            Cqu = zeros(P.Model.NDof,P.Mesh.NExcite,NPts);
            Cxu = zeros(P.Model.NDofInt,P.Mesh.NExcite,NPts);
            u0 = u;
            udot0 = udot;
            for i = 1:P.Mesh.NExcite
                u(i,:) = u(i,:) + h;
                f = rotor_hbm_model('all',t,x,xdot,xddot,u,udot,uddot,hbm,problem,w0);
                Kqu(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Kxu(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                u = u0;
                
                udot(i,:) = udot(i,:) + h;
                f = rotor_hbm_model('all',t,x,xdot,xddot,u,udot,uddot,hbm,problem,w0);
                Cqu(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Cxu(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                udot = udot0;
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