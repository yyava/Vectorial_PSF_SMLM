function [logL,gradlogL,HessianlogL] = likelihood(p,image,mu,dmu)
% Returns the log-likelihood, as well as first and second order
% derivatives w.r.t. the parameters for a noisy image M, measured in
% number of photons per pixel, and Poisson-rate mu with first order
% derivatives dmudtheta. The log-likelihood is Poisson+readout-noise based.
% Input: 
%   image(Nx,Ny,Nz,Nc)
%   mu(Nx,Ny,Nz,Nc)
%   dmu_dTheta(Nx,Ny,Nz,Nc,Np)
% Output: 
%   logL(1)
%   gradlogL(Np,1)
%   HessianlogL(Np,Np)

Np = p.Np;
varfit = p.varfit;    % Gaussian readout noise

% calculation of weight factors
keps = 1e3*eps;
mupos = double(mu>0).*mu + double(mu<=0)*keps;

weight = (image-mupos)./(mupos+varfit);
dweight = (image+varfit)./(mupos+varfit).^2;

% log-likelihood, gradient vector and Hessian matrix
logL = sum((image+varfit).*log(mupos+varfit)-(mupos+varfit),"all");

gradlogL = reshape(sum(weight.*dmu,1:4),[Np,1]);

HessianlogL = zeros(Np,Np);
for ii = 1:Np
    for jj = ii:Np
        HessianlogL(ii,jj) = sum(-dweight.*dmu(:,:,:,:,ii).*dmu(:,:,:,:,jj),1:4);
        HessianlogL(jj,ii) = HessianlogL(ii,jj);
    end
end




