function [mu,dmu] = get_PoissonRate(p, Theta)
% This function calculates the excepted photon counts (Poisson-rates) and 
% its derivative with given fitting parameter Theta
% Input: 
%   Theta(Np)
% Output: 
%   mu(Nx,Ny,1,Nc)
%   dmu_dTheta(Nx,Ny,1,Nc,Np) Np=5-8
%   Np:[dx,dy,dz,d_Nph,d_Nbg,dazim,dpola,dg2]

fitModel=p.fitModel;

p.x0 = Theta(1);
p.y0 = Theta(2);
p.z0 = Theta(3);
p.Nph = Theta(4);
p.Nbg = Theta(5);

if contains(fitModel,'azim-pola')
    p.azim = Theta(6);
    p.pola = Theta(7);
end
if contains(fitModel,'diffusion')
    p.g2 = Theta(8);
end

[PupilMatrix,dPupilMatrix,Energy_norm] = get_pupilMatrix(p);
[PSF,dPSF] = get_PSF(p,PupilMatrix,dPupilMatrix,Energy_norm);
[mu,dmu] = get_mu(p,PSF,dPSF); 
end