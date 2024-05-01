function [mu,dmu] = get_mu(p, PSF, dPSF)
% This function calculates the excepted photon counts (Poisson-rates) and 
% its derivative
% Input: 
%   PSF size:(Nx,Ny,Nz,Nc)
%   dPSF size:(Nx,Ny,Nz,Nc,Np) Np=3-6
%       D5(dPSF): derivatives [dx,dy,dz(,dazim,dpola,dg2)]
% Output: 
%   mu size:(Nx,Ny,Nz,Nc)
%   dmu size:(Nx,Ny,Nz,Nc,Np) Np=5-8
%       D5(dmu): derivatives [dx,dy,dz,d_Nph,d_Nbg(,dazim,dpola,dg2)]

% signal and background photon counts
Nph=p.Nph;
Nbg=p.Nbg/p.Nc; % mean Nbg in all channels
fitModel=p.fitModel;

if ~p.Excitation
    % excepted photon counts
    if ~strcmp(p.polarization,'None')
        mu = Nph*PSF+Nbg;
    else
        mu = Nph*PSF;
        mu(:,:,:,1) = mu(:,:,:,1)+p.Nbg;
    end
    % Derivatives
    % [dx,dy,dz]
    dmu(:,:,:,:,1:3) = Nph*dPSF(:,:,:,:,1:3);
    % [d_Nph]
    dmu(:,:,:,:,4) = PSF;
    % [d_Nbg]
    if ~strcmp(p.polarization,'None')
        dmu(:,:,:,:,5) = 1/p.Nc;
    else
        dmu(:,:,:,1,5) = 1;
        dmu(:,:,:,2:4,5) =0;
    end
    if contains(fitModel,'azim-pola')
        % [dazim,dpola]
        dmu(:,:,:,:,6:7) = Nph*dPSF(:,:,:,:,4:5);
    end
    if contains(fitModel,'diffusion')
        % [dg2]
        dmu(:,:,:,:,8) = Nph*dPSF(:,:,:,:,6);
    end
else    % Excitation modulation
    M = p.M;
    Xi = p.Xi;
    mu = zeros(p.Nx,p.Nx,p.Nz,p.Nc*M);
    dmu = zeros(p.Nx,p.Nx,p.Nz,p.Nc*M,p.Np);
    for j = 1:M
        NphNorm = Nph*2*cos(p.azim-Xi(j)).^2/M;
        mu(:,:,:,(j-1)*p.Nc+1:j*p.Nc) = (NphNorm*PSF+Nbg/M);
        dmu(:,:,:,(j-1)*p.Nc+1:j*p.Nc,1:3) = NphNorm*dPSF(:,:,:,:,1:3);
        % [dNph]
        dmu(:,:,:,(j-1)*p.Nc+1:j*p.Nc,4) = NphNorm*PSF/Nph;
        if contains(fitModel,'azim-pola')
        % [dazim,dpola]
            dmu(:,:,:,(j-1)*p.Nc+1:j*p.Nc,6) = ...
            (dPSF(:,:,:,:,4)*2*cos(p.azim-Xi(j)).^2-4*PSF*cos(p.azim-Xi(j))*sin(p.azim-Xi(j)))*Nph/M;
            dmu(:,:,:,(j-1)*p.Nc+1:j*p.Nc,7) = NphNorm*dPSF(:,:,:,:,5);
        end
        if contains(fitModel,'diffusion')
            % [dg2]
            dmu(:,:,:,(j-1)*p.Nc+1:j*p.Nc,8) = NphNorm*dPSF(:,:,:,:,6);
        end
    end
    % [dNbg]
    if ~strcmp(p.polarization,'None')
        dmu(:,:,:,:,5) = 1/p.Nc/M;
    else
        dmu(:,:,:,:,5) = 0;
        dmu(:,:,:,[1,1+p.Nc,1+2*p.Nc],5) = 1/M;
    end
end


