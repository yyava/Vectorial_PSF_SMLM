function theta_init = initValues(allspots,p)
% Provide initial values for parameter fitting
% Input: 
%   allspots(Nx,Ny,Nz,Nc,Ncfg)
% Output: 
%   theta_init(Np,Ncfg)
%   Np:[x0,y0,z0,Nph,Nbg,azim,pola,g2]

Nz=p.Nz; Nx=p.Nx; Np=p.Np; Nc=p.Nc; Ncfg=p.Ncfg;
fitModel=p.fitModel;

% sampling coordinates in image plane
[XX,YY] = meshgrid(p.xl,p.xl);
photonFlux = 1.5; % rough photon flux correction

theta_init = zeros(Np,Ncfg);

azim0 = pi/4;
pola0 = pi/3;
g20 = 0.75;

for jcfg = 1:Ncfg
    % rename config
    fxyz = sum(allspots(:,:,:,:,jcfg),4);
    fxy = fxyz(:,:,ceil(Nz/2));
    % background estimate from the median value of the rim pixels
    fxytmp = fxy';
    rimpx = [fxytmp(1:Nx) fxytmp(end-(Nx-1):end) fxy(1:Nx) fxy(end-(Nx-1):end)];
    bg = min(max(median(rimpx),1),50);
    % estimate of signal photon count
    fxy = fxy-bg;
    % raw moments
    m00 = sum(fxy,1:2);
    m10 = sum(XX.*fxy,1:2);
    m01 = sum(YY.*fxy,1:2);
    % m20 = sum(XX.^2.*fxy,1:2);
    % m02 = sum(YY.^2.*fxy,1:2);
    % m11 = sum(XX.*YY.*fxy,1:2);
    % centroids (lateral)
    xc = m10/m00;
    yc = m01/m00;
    
    % estimate axial position
    z0 = 0e-9;

    theta_init(1:5,jcfg) = [xc,yc,z0,m00*photonFlux,bg/Nc];
    if contains(fitModel,'azim-pola')
        theta_init(6:7,jcfg) = [azim0,pola0];
    end
    if contains(fitModel,'diffusion')
        theta_init(8,jcfg) = g20;
    end

end

end