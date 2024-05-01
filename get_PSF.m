function [PSF,dPSF] = get_PSF(p,PupilMatrix,dPupilMatrix)
% Returns to the normalized PSF of the dipole emitter as the
% superposition of both fixed and free dipole.
% Input: 
%   PupilMatrix size:(Nk,Nk,Nz,Nc,3)
%   dPupilMatrix size:(Nk,Nk,Nz,Nc,3,3)
%       D6(dPupilMatrix): spatial derivatives [dx,dy,dz]
% Output: 
%   PSF size:(Nx,Ny,Nz,Nc)
%   dPSF size:(Nx,Ny,Nz,Nc,Np) Np=3-6
%       D5(dPSF): derivatives [dx,dy,dz(,dazim,dpola,dg2)]

Nx=p.Nx; Nz=p.Nz; Nc=p.Nc; 
pola=p.pola; azim=p.azim;g2=p.g2;
fitModel=p.fitModel;
% czt parameters
M = p.Nx; w = p.czt_w; a = p.czt_a;

% dipole orientation unit vector
dipoVec=[sin(pola)*cos(azim),sin(pola)*sin(azim),cos(pola)];

if strcmp(p.dipoleType,'fixed')
    % FT after multiplying of electric field and dipole vector
    [PSF,dPSF] = get_PSF_fixed_4pi(p,PupilMatrix,dPupilMatrix,dipoVec);
else
    % electric field and its derivatives of 3 vector components in 4 channels
    E = czt2D(PupilMatrix,M,w,a)*(p.Dx*p.Dk/2/pi);    % size:(Nx,Ny,Nz,Nc,3(xyz))
    dEdr = czt2D(dPupilMatrix,M,w,a)*(p.Dx*p.Dk/2/pi);  % size:(Nx,Ny,Nz,Nc,3(xyz),3(dr))
    % PSF and its derivative for free dipole
    PSF_free = 1/3*sum(abs(E).^2,5);
    dPSFdr_free = 2/3*real(reshape(sum(dEdr.*conj(E),5),[Nx,Nx,Nz,Nc,3]));
    % normalization, finding the total energy in pupil plane
    Energy_free = 1/3*sum(abs(PupilMatrix).^2,[1,2,4,5]);
    Energy_fixed = sum(abs(...
        PupilMatrix(:,:,:,:,1)*dipoVec(1)+...
        PupilMatrix(:,:,:,:,2)*dipoVec(2)+...
        PupilMatrix(:,:,:,:,3)*dipoVec(3)).^2,[1,2,4]);
if strcmp(p.dipoleType,'free')
    PSF = PSF_free./Energy_free;
    dPSF = dPSFdr_free./Energy_free;
elseif strcmp(p.dipoleType,'diffusion')
    % electric field of fixed dipole on image plane
    E_fixed = E(:,:,:,:,1)*dipoVec(1)+...
        E(:,:,:,:,2)*dipoVec(2)+E(:,:,:,:,3)*dipoVec(3);
    % PSF of fixed dipole
    PSF_fixed = abs(E_fixed).^2;
    % total PSF as the composition of fixed and free dipole
    PSF = (1-g2)*PSF_free+g2*PSF_fixed;
    dEdr_fixed = reshape((dEdr(:,:,:,:,1,:)*dipoVec(1)+dEdr(:,:,:,:,2,:)...
    *dipoVec(2)+dEdr(:,:,:,:,3,:)*dipoVec(3)),[Nx,Nx,Nz,Nc,3]);
    dPSFdr_fixed = 2*real(dEdr_fixed.*conj(E_fixed));
    % spatial derivative of PSF
    dPSF = (1-g2)*dPSFdr_free+g2*dPSFdr_fixed;
    if contains(fitModel,'azim-pola')     
        % angular derivatives of dipole vector 
        dazim = [-sin(pola)*sin(azim),sin(pola)*cos(azim),0];
        dpola = [cos(pola)*cos(azim),cos(pola)*sin(azim),-sin(pola)];
        dEdazim_fixed = E(:,:,:,:,1)*dazim(1)+E(:,:,:,:,2)*dazim(2)+E(:,:,:,:,3)*dazim(3);
        dEdpola_fixed = E(:,:,:,:,1)*dpola(1)+E(:,:,:,:,2)*dpola(2)+E(:,:,:,:,3)*dpola(3);
        % [dazim,dpola]
        dPSF(:,:,:,:,4) = g2*2*real(dEdazim_fixed.*conj(E_fixed));
        dPSF(:,:,:,:,5) = g2*2*real(dEdpola_fixed.*conj(E_fixed));
    end
    if contains(fitModel,'diffusion')
    % [dg2]
    dPSF(:,:,:,:,6) = PSF_fixed-PSF_free;
    end
    % Normalization
    tot_E = (1-p.g2)*Energy_free+p.g2*Energy_fixed;
    PSF = PSF./tot_E;
    dPSF = dPSF./tot_E;
end
end

% Non-polarized
if strcmp(p.polarization,'None')
    PSF(:,:,:,1) = sum(PSF,4);
    PSF(:,:,:,2:4) = 0;
    dPSF(:,:,:,1,:) = sum(dPSF,4);
    dPSF(:,:,:,2:4,:) = 0;    
end

end

%% local functions
function [PSF_fixed,dPSF_fixed] = get_PSF_fixed_4pi(p,PupilMatrix,dPupilMatrix,dipoVec)
% This function calculates the PSF and its derivative for fixed dipole
% case. This function skips calculating the electric field of 3 vector
% components, so only a third of the number of FTs are needed
    pola=p.pola; azim=p.azim;
    Nz=p.Nz; Nc=p.Nc; Nk=p.Nk;
    fitModel=p.fitModel;
    M = p.Nx; w = p.czt_w; a = p.czt_a;
    % 
    PupilMatrix_fixed = PupilMatrix(:,:,:,:,1)*dipoVec(1)...
        +PupilMatrix(:,:,:,:,2)*dipoVec(2)+PupilMatrix(:,:,:,:,3)*dipoVec(3);
    E_fixed = czt2D(PupilMatrix_fixed,M,w,a)*(p.Dx*p.Dk/2/pi);
    PSF_fixed = abs(E_fixed).^2;
    dPupilMatrixdr_fixed = reshape(dPupilMatrix(:,:,:,:,1,:)*dipoVec(1)...
        +dPupilMatrix(:,:,:,:,2,:)*dipoVec(2)...
        +dPupilMatrix(:,:,:,:,3,:)*dipoVec(3),[Nk,Nk,Nz,Nc,3]);
    dEdr_fixed = czt2D(dPupilMatrixdr_fixed,M,w,a)*(p.Dx*p.Dk/2/pi);
    dPSF_fixed = 2*real(dEdr_fixed.*conj(E_fixed));
    if contains(fitModel,'azim-pola')
        % angular derivatives of dipole vector 
        dazim = [-sin(pola)*sin(azim),sin(pola)*cos(azim),0];
        dpola = [cos(pola)*cos(azim),cos(pola)*sin(azim),-sin(pola)];
        dPupilMatrixdazim_fixed = PupilMatrix(:,:,:,:,1)*dazim(1)...
            +PupilMatrix(:,:,:,:,2)*dazim(2)+PupilMatrix(:,:,:,:,3)*dazim(3);
        dPupilMatrixdpola_fixed = PupilMatrix(:,:,:,:,1)*dpola(1)...
            +PupilMatrix(:,:,:,:,2)*dpola(2)+PupilMatrix(:,:,:,:,3)*dpola(3);
        dEdazim_fixed = czt2D(dPupilMatrixdazim_fixed,M,w,a)*(p.Dx*p.Dk/2/pi);
        dEdpola_fixed = czt2D(dPupilMatrixdpola_fixed,M,w,a)*(p.Dx*p.Dk/2/pi);
        dPSF_fixed(:,:,:,:,4) = 2*real(dEdazim_fixed.*conj(E_fixed));
        dPSF_fixed(:,:,:,:,5) = 2*real(dEdpola_fixed.*conj(E_fixed));
    end
    Energy_fixed = sum(abs(PupilMatrix_fixed).^2,[1,2,4]);
    PSF_fixed = PSF_fixed./Energy_fixed;
    dPSF_fixed = dPSF_fixed./Energy_fixed;
end
