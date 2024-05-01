function [CRLB,FisherM] = get_CRLB(p,mu,dmu)
% Returns to the CRLB and Fisher Matrix of all parameters
% Input: 
%   mu(Nx,Ny,Nz,Nc)
%   dmu(Nx,Ny,Nz,Nc,Np) Np=5-8
% Output: 
% zslice:
%   CRLB(Np)
%   FisheM(Np,Np)
% zstack:
%   CRLB(Np,Nz)
%   FisheM(Np,Np,Nz)
% Np:[dx,dy,dz,d_Nph,d_Nbg,dazim,dpola,dg2]

keps = 1e3*eps;
warning('off','all')


%% parameter setting
Np = p.Np;
Nz = p.Nz;

mupos = (mu>0).*mu+(mu<=0).*keps;   % Avoid dividing by zero
if strcmp(p.detection,"zslice")
    FisherM = zeros(Np,Np);
    for ii = 1:Np
    for jj = ii:Np
        FisherM(ii,jj) = sum(1./mupos.*dmu(:,:,:,:,ii).*dmu(:,:,:,:,jj),"all");
        FisherM(jj,ii) = FisherM(ii,jj);
    end
    end
    CRLB = sqrt(diag(inv(FisherM(:,:)+keps*eye(Np))));

elseif strcmp(p.detection,"zstack")
    FisherM = zeros(Np,Np,Nz);
    for ii = 1:Np
        for jj = ii:Np
            FisherM(ii,jj,:) = sum(1./mupos.*dmu(:,:,:,:,ii).*dmu(:,:,:,:,jj),[1,2,4]);
            FisherM(jj,ii,:) = FisherM(ii,jj,:);
        end
    end
    CRLB = zeros(Np,Nz);
    for iz = 1:Nz
        CRLB(:,iz) = sqrt(diag(inv(FisherM(:,:,iz)+keps*eye(Np))));
    end
end
