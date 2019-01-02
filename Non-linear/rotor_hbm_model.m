function varargout = rotor_hbm_model(part,t,x,xdot,xddot,u,udot,uddot,xalg,hbm,problem,w0)
P = problem.P;
O  = w0(1)/hbm.harm.rFreqRatio(1);
NPts = size(x,2);
tRot = t(P.Model.iRot,:);

NDof = size(x,1);
NInput = size(u,1);

if P.Model.bUseAlgebraic
    xInt = xalg;
    
    xCG = x;
    xdotCG = xdot;
    xddotCG = xddot;
else
    xInt = x(P.Model.NDof+1:end,:);

    xCG = x(1:P.Model.NDof,:);
    xdotCG = xdot(1:P.Model.NDof,:);
    xddotCG = xddot(1:P.Model.NDof,:);
end

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

if P.Model.bUseAlgebraic
    xalg = States.xInt;
else
    x = [xCG; States.xInt];
    xdot = [xdotCG; States.xdotInt];
    xddot = [xddotCG; States.xddotInt];
end

varargout = {};
switch part
    case  {'nl','alg','all','output'}
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
        elseif strcmp(part,'alg')
            if P.Model.bUseAlgebraic
                varargout{end+1} = Fi;
            else
                varargout{end+1} = zeros(0,NPts);
            end
        elseif strcmp(part,'nl')
            if P.Model.bUseAlgebraic
                varargout{end+1} = Fnl;
            else
                varargout{end+1} = [Fnl; Fi];
            end
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
            f0 = rotor_hbm_model('all',t,x,xdot,xddot,u,udot,uddot,xalg,hbm,problem,w0);
            h = 1E-10;
            
            Kqq = zeros(P.Model.NDof,P.Model.NDof,NPts);
            Kxq = zeros(P.Model.NDofInt,P.Model.NDof,NPts);
            Cqq = zeros(P.Model.NDof,P.Model.NDof,NPts);
            Cxq = zeros(P.Model.NDofInt,P.Model.NDof,NPts);
            x0 = x;
            xdot0 = xdot;
            for i = 1:P.Model.NDof
                x(i,:) = x(i,:) + h;
                f = rotor_hbm_model('all',t,x,xdot,xddot,u,udot,uddot,xalg,hbm,problem,w0);
                Kqq(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Kxq(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                x = x0;
                
                xdot(i,:) = xdot(i,:) + h;
                f = rotor_hbm_model('all',t,x,xdot,xddot,u,udot,uddot,xalg,hbm,problem,w0);
                Cqq(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Cxq(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                xdot = xdot0;
            end
            
            Kqx = zeros(P.Model.NDof,P.Model.NDofInt,NPts);
            Kxx = zeros(P.Model.NDofInt,P.Model.NDofInt,NPts);
            Cqx = zeros(P.Model.NDof,P.Model.NDofInt,NPts);
            Cxx = zeros(P.Model.NDofInt,P.Model.NDofInt,NPts);
            xalg0 = xalg;
            for i = 1:P.Model.NDofInt
                if P.Model.bUseAlgebraic
                    xalg(i,:) = xalg(i,:) + h;
                else
                    x(i + P.Model.NDof,:) = x(i + P.Model.NDof,:) + h;
                end
                f = rotor_hbm_model('all',t,x,xdot,xddot,u,udot,uddot,xalg,hbm,problem,w0);
                Kqx(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Kxx(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                x = x0;
                xalg = xalg0;
                
                Cqx(:,i,:) = zeros(P.Model.NDof,1,NPts);
                Cxx(:,i,:) = zeros(P.Model.NDofInt,1,NPts);
            end
            
            Kqu = zeros(P.Model.NDof,P.Mesh.NExcite,NPts);
            Kxu = zeros(P.Model.NDofInt,P.Mesh.NExcite,NPts);
            Cqu = zeros(P.Model.NDof,P.Mesh.NExcite,NPts);
            Cxu = zeros(P.Model.NDofInt,P.Mesh.NExcite,NPts);
            u0 = u;
            udot0 = udot;
            for i = 1:P.Mesh.NExcite
                u(i,:) = u(i,:) + h;
                f = rotor_hbm_model('all',t,x,xdot,xddot,u,udot,uddot,xalg,hbm,problem,w0);
                Kqu(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Kxu(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                u = u0;
                
                udot(i,:) = udot(i,:) + h;
                f = rotor_hbm_model('all',t,x,xdot,xddot,u,udot,uddot,xalg,hbm,problem,w0);
                Cqu(:,i,:) = (f(1:P.Model.NDof,:)-f0(1:P.Model.NDof,:))/h;
                Cxu(:,i,:) = (f(P.Model.NDof+1:end,:)-f0(P.Model.NDof+1:end,:))/h;
                udot = udot0;
            end
        end
        
        switch part
            case {'nl_x','nl_alg','alg_x','alg_alg'}
                if bReduced
                    [Kqq,Kqx,Kxq,Kxx] = compress_state_derivatives(P,O,hbm,Kqq,Kqx,Kxq,Kxx);
                end
     
                if P.Model.bUseAlgebraic
                    switch part
                        case 'nl_x'
                            varargout{end+1} = Kqq;
                        case 'nl_alg'
                            varargout{end+1} = Kqx;
                        case 'alg_x'
                            varargout{end+1} = Kxq;
                        case 'alg_alg'
                            varargout{end+1} = Kxx;
                    end
                else
                    switch part
                        case 'nl_x'
                            varargout{end+1} = [Kqq Kqx; Kxq Kxx];
                        case 'nl_alg'
                            varargout{end+1} = zeros(NDof,0,NPts);
                        case 'alg_x'
                            varargout{end+1} = zeros(0,NDof,NPts);
                        case 'alg_alg'
                            varargout{end+1} = zeros(0,0,NPts);
                    end
                end
            case {'nl_xdot','alg_xdot'} 
                if bReduced
                    [Cqq,Cqx,Cxq,Cxx] = compress_state_derivatives(P,O,hbm,Cqq,Cqx,Cxq,Cxx);
                end
                if P.Model.bUseAlgebraic
                    switch part
                        case 'nl_xdot'
                            varargout{end+1} = Cqq;
                        case 'alg_xdot'
                            varargout{end+1} = Cxq;
                    end
                else
                    switch part
                        case 'nl_xdot'
                            varargout{end+1} = [Cqq Cqx; Cxq Cxx];
                        case 'alg_xdot'
                            varargout{end+1} = zeros(0,NDof,NPts);
                    end
                end
            case {'nl_u','alg_u'}
                if bReduced
                    [Kqu,Kxu] = compress_input_derivatives(P,O,hbm,Kqu,Kxu);
                end
                 if P.Model.bUseAlgebraic
                     switch part
                         case 'nl_u'
                             varargout{end+1} = Kqu;
                         case 'alg_u'
                             varargout{end+1} = Kxu;
                     end
                else
                    switch part
                        case 'nl_u'
                            varargout{end+1} = [Kqu; Kxu];
                        case 'alg_u'
                            varargout{end+1} = zeros(0,NInput,NPts);
                    end
                 end
            case {'nl_udot','alg_udot'}
                if bReduced
                    [Cqu,Cxu] = compress_input_derivatives(P,O,hbm,Cqu,Cxu);
                end
                if P.Model.bUseAlgebraic
                     switch part
                         case 'nl_udot'
                             varargout{end+1} = Cqu;
                         case 'alg_udot'
                             varargout{end+1} = Cxu;
                     end
                else
                    switch part
                        case 'nl_udot'
                            varargout{end+1} = [Cqu; Cxu];
                        case 'alg_udot'
                            varargout{end+1} = zeros(0,NInput,NPts);
                    end
                end
            case 'nl_xddot'
                if P.Model.bUseAlgebraic
                    varargout{end+1} = 0*Kqx;
                else
                    varargout{end+1} = 0*[Kqx; Kxx];
                end
            case 'alg_xddot'
                if P.Model.bUseAlgebraic
                    varargout{end+1} = 0*Kxx;
                else
                    varargout{end+1} = zeros(0,NDof,NPts);
                end
            case 'nl_uddot'
                if P.Model.bUseAlgebraic
                     varargout{end+1} = 0*Kqu;
                else
                    varargout{end+1} = 0*[Kqu; Kxu];
                end
            case 'alg_uddot'
                if P.Model.bUseAlgebraic
                    varargout{end+1} = 0*Kxu;
                else
                    varargout{end+1} = zeros(0,NInput,NPts);
                end
            case {'nl_w' 'alg_w'}
                varargout{end+1} = [];
                %too complicated to do analytically
        end
end

function xInt = saturate_bearing_states(P,O,hbm,xInt)
iDofOut = 0;
for i = 1:length(P.Bearing)
    for j = 1:2
        if strcmp(P.Bearing{i}.Model{j},'REB')
            REB = P.Bearing{i}.Params{j};

            Ocage = REB.Kinematics.rCagei*O;
            
            Fc = (REB.Dynamics.mb * Ocage.^2 * REB.Geometry.dm/2);
            xc = (Fc/REB.Contact.Outer.K).^(1/REB.Contact.Outer.n);
            xInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:) = max(xInt(iDofOut+(1:P.Bearing{i}.NDofInt(j)),:),xc);
            iDofOut = iDofOut + P.Bearing{i}.NDofInt(j);
        else
            
            
        end
    end
end

function [Kqq2,Kqx2,Kxq2,Kxx2] = compress_state_derivatives(P,O,hbm,Kqq,Kqx,Kxq,Kxx)
Red = P.Model.Reduced;

Kqq2 = Kqq;
Kqx2 = Kqx(:,Red.iInt,:);
Kxq2 = Kxq(Red.iInt,:,:);
Kxx2 = Kxx(Red.iInt,Red.iInt,:);

function [Kqu2,Kxu2] = compress_input_derivatives(P,O,hbm,Kqu,Kxu)
Red = P.Model.Reduced;
Kqu2 = Kqu;
Kxu2 = Kxu(Red.iInt,:,:);