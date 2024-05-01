function [ThetaStore,muStore,dmuStore,meritStore,numiters,ThetaAllStore] = localization(allspots,Theta0,p)
fprintf('\nStart fitting %i instances:\n',p.Ncfg); tic;
% This function finds the parameters for a single 2D-image assuming
% vectorial PSF models.
% Output: 
%   ThetaStore(Np,NiterMax+1,Ncfg)
%   muStore(Nx,Nx,Nz,Nc,Ncfg)
%   dmuStore(Nx,Nx,Nz,Nc,Np,Ncfg)
%   meritStore(NiterMax+1,Ncfg)
%   numiters(Ncfg,1)

% parameter settings
Ncfg = p.Ncfg;
tollim = p.tollim;
NiterMax = p.NiterMax;
Np=p.Np;
Nx=p.Nx;
Nz=p.Nz;
Nc=p.Nc;

% pre-allocation
numiters = zeros(Ncfg,1);
muStore = zeros(Nx,Nx,Nz,Nc,Ncfg);
dmuStore = zeros(Nx,Nx,Nz,Nc,Np,Ncfg);
ThetaAllStore = zeros(Np,NiterMax+1,Ncfg);
meritStore = zeros(NiterMax+1,Ncfg);
if p.Excitation
    muStore = zeros(Nx,Nx,Nz,Nc*p.M,Ncfg);
    dmuStore = zeros(Nx,Nx,Nz,Nc*p.M,Np,Ncfg);
end

% setup parallel loop
if p.flg_parallel
    if isempty(gcp('nocreate'))
        parpool;
    end
    parforArg = Inf;
    fprintf('parallel computation start...\n')
else
    parforArg = 0;
end

parfor (jcfg = 1:Ncfg, parforArg)
    % pre-allocate
    ThetaTemp = ones(Np,NiterMax+1);
    meritTemp = zeros(NiterMax+1,1);
    
    % initial values  
    Theta = Theta0(:,jcfg);
    spots = allspots(:,:,:,:,jcfg);

    [mu,dmu] = get_PoissonRate(p,Theta);
    [merit,grad,Hessian] = likelihood(p,spots,mu,dmu);

    ThetaTemp(:,1) = Theta;
    meritTemp(1) = merit;

    % start iteration loop
    iiter = 1; 
    monitor = 2*tollim;
    alambda = 1;
    alambdafac = 10;
    
    while ((iiter<=NiterMax) && (monitor>tollim))
        % update parameters
        ThetaNew = ThetaUpdate(Theta,grad,Hessian,alambda,p);
        % calculate update merit function
        [mu,dmu] = get_PoissonRate(p,ThetaNew);
        [meritTry,gradTry,HessianTry] = likelihood(p,spots,mu,dmu);

        % modify Levenberg-Marquardt parameter
        if (meritTry>merit)
            monitor = abs((meritTry-merit)/merit);
            alambda = alambda/alambdafac;
            Theta = ThetaNew;
            merit = meritTry;
            grad = gradTry;
            Hessian = HessianTry;
        else
            alambda = alambdafac*alambda;
        end
        ThetaTemp(:,iiter+1) = Theta;
        meritTemp(iiter+1) = merit;
        iiter = iiter+1;
    end
    
    % store values
    meritTemp(iiter+1:end) = merit;
    ThetaTemp(:,iiter+1:end) = ThetaTemp(:,iiter+1:end).*Theta;

    ThetaAllStore(:,:,jcfg) = ThetaTemp;
    muStore(:,:,:,:,jcfg) = mu;
    dmuStore(:,:,:,:,:,jcfg) = dmu;
    meritStore(:,jcfg) = meritTemp;
    numiters(jcfg) = iiter;

    % print update
    if rem(jcfg,round(Ncfg/10)) == 0
        fprintf('fitting instance # %i...\n',jcfg)
    end

    % add offset to merit function for an accurate determination of the log-likelihood
    % meritOffset = -sum(gammaln(allspots(:,:,:,:,jcfg)+p.varfit+1),"all");
    % meritStore(:,jcfg) = meritStore(:,jcfg)+meritOffset;
end

ThetaStore = squeeze(ThetaAllStore(:,end,:));

% print run time
fprintf(['\nMLE fit routine (spot/second): ' num2str(toc,3) 's (' num2str(p.Ncfg/toc,5) ')\n'])
