function ThetaNew = ThetaUpdate(ThetaOld,grad,Hessian,alambda,p)
% This function calculates the new parameters in each iteration step.
%
warning('off','all')

fitModel = p.fitModel;

% update of fit parameters via Levenberg-Marquardt
Bmat = Hessian+alambda*diag(diag(Hessian));
if (abs(det(Bmat))<1e3*eps)
    Bmat = Bmat+1e3*eps*eye(size(Bmat));
end

dTheta = -Bmat\grad;

ThetaNew = ThetaOld+dTheta;

% enforce physical boundaries in angular space
if contains(fitModel,'azim-pola')
    ThetaNew(6) = mod(ThetaNew(6),2*pi);
    ThetaNew(7) = mod(ThetaNew(7),pi);
    if ThetaNew(7)>pi/2
        ThetaNew(7) = pi-ThetaNew(7);
        ThetaNew(6) = mod(ThetaNew(6)+pi,2*pi);
    end
end
if contains(fitModel,'diffusion')
    if ThetaNew(8)-ThetaOld(8)<-0.1
        ThetaNew(8) = ThetaOld(8)-0.1;
    end
    if ThetaNew(8)>1.0
        ThetaNew(8) = 1.0;
    elseif ThetaNew(8)<0.0
        ThetaNew(8) = 0.0;
    end
end

% enforce physical boundaries in parameter space.
if abs(ThetaNew(1))>p.Lx/4
    ThetaNew(1) = ThetaOld(1);
end
if abs(ThetaNew(2))>p.Lx/4
    ThetaNew(2) = ThetaOld(2);
end
if ThetaNew(3)<p.zrange(1) || ThetaNew(3)>p.zrange(2)
    ThetaNew(3) = ThetaOld(3);
    % ThetaNew(3) = 0;
end

