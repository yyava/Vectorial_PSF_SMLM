function [PupilMatrix,dPupilMatrix,Energy_norm] = get_pupilMatrix(p)
% This fuction calculates the pupil matrix its derivative
% Output: 
%   PupilMatrix size:(Nk,Nk,Nz,Nc,3)
%   dPupilMatrix size:(Nk,Nk,Nz,Nc,3,3)
%       D1~3: kx,ky,z
%       D4: 4 channels [x1,x2,y1,y2]
%       D5: vector components [x,y,z]
%       D6(dPupilMatrix): spatial derivatives of x,y,z [dx,dy,dz]

Nk=p.Nk; Nz=p.Nz; Nc=p.Nc;  
pupila=p.pupila; pupilb=p.pupilb;   % aberrations included 
cp = p.chPhase;     % phase shift of arm b in 4 channels
zl = p.zl+p.z0;     % axial shift
if p.NA>p.n_med     % SAF
    zl = zl+p.depth;     % axial image plane position with shift
end

pupilMask = p.pupilMask;
k0 = p.k0; kr = p.kr; phi = p.phi;
n_med = p.n_med; n_cov = p.n_cov; n_imm = p.n_imm;
kx=p.kx; ky=p.ky;   % wavevector components
kz=p.kz; kzimm=p.kzimm;
x0=p.x0; y0=p.y0;   % lateral emitter position

% Wavevector z-component of the incident medium, this factor originates from 
% the Weyl-representation of the emitted vector spherical wave of the dipole.
% theta in different media // sin_med means sin(\theta_med)
sin_med = kr/k0/n_med.*pupilMask;
sin_cov = kr/k0/n_cov.*pupilMask;
sin_imm = kr/k0/n_imm.*pupilMask;
cos_med = sqrt(1-sin_med.^2);
cos_cov = sqrt(1-sin_cov.^2);
cos_imm = sqrt(1-sin_imm.^2);

% Fresnel coefficients for p and s polarization
% transmission of interfaces: medium-coverslip-immersion fluid
Tp = 4*n_med*n_cov*cos_med.*cos_cov./...
    (n_cov*cos_med+n_med*cos_cov)./...
    (n_cov*cos_imm+n_imm*cos_cov);
Ts = 4*n_med*n_cov*cos_med.*cos_cov./...
    (n_med*cos_med+n_cov*cos_cov)./...
    (n_cov*cos_cov+n_imm*cos_imm);

% polarization vectors
p_unit(:,:,1) = cos_med.*cos(phi).*pupilMask;
p_unit(:,:,2) = cos_med.*sin(phi).*pupilMask;
p_unit(:,:,3) = -sin_med.*pupilMask;
s_unit(:,:,1) = -sin(phi).*pupilMask;
s_unit(:,:,2) = cos(phi).*pupilMask;
s_unit(:,:,3) = zeros(p.Nk).*pupilMask;
p_vec = p_unit.*Tp;
s_vec = s_unit.*Ts;

% polarization vectors along x y axis
qX = cos(phi).*p_vec-sin(phi).*s_vec;
qY = sin(phi).*p_vec+cos(phi).*s_vec;
% 4LP
qXr = (qX+qY)/sqrt(2);
qYr = (qX-qY)/sqrt(2);
% Circular polarization
qL = (qX+1i*qY)/sqrt(2);
qR = (1i*qX+qY)/sqrt(2);

% Elliptical polarization
chi = p.chi; % rotation angle of QWP (elipticity)

Rp = [cosd(chi),sind(chi);-sind(chi),cosd(chi)];
Rn = [cosd(-chi),sind(-chi);-sind(-chi),cosd(-chi)];

QWPp = Rn*[1,0;0,1i]*Rp;
QWPn = Rp*[1,0;0,1i]*Rn;

qXe = QWPp(1)*qX+QWPp(2)*qY;
qYe = QWPp(3)*qX+QWPp(4)*qY;
% For chi=22.5°
qXre = QWPn(1)*qX+QWPn(2)*qY;
qYre = QWPn(3)*qX+QWPn(4)*qY;

% Arbitary polarization using composite wave-plate (CEP)
% HWP = sqrt(0.5)*[1,1;1,-1];
% CWP = HWP*QWPn;
% qXre = CWP(1)*qX+CWP(2)*qY;
% qYre = CWP(3)*qX+CWP(4)*qY;

% matching variable index for 4 channels
if strcmp(p.polarization,'2LP')||strcmp(p.polarization,'None')
    PolaVecA(:,:,:,1) = qX;    PolaVecA(:,:,:,2) = qY;
    PolaVecA(:,:,:,3) = qX;    PolaVecA(:,:,:,4) = qY;
elseif strcmp(p.polarization,'4LP')
    PolaVecA(:,:,:,1) = qX;    PolaVecA(:,:,:,2) = qY;
    PolaVecA(:,:,:,3) = qXr;   PolaVecA(:,:,:,4) = qYr;
elseif strcmp(p.polarization,'CP')
    PolaVecA(:,:,:,1) = qL;    PolaVecA(:,:,:,2) = qR;
    PolaVecA(:,:,:,3) = qL;    PolaVecA(:,:,:,4) = qR;
elseif strcmp(p.polarization,'2EP')
    PolaVecA(:,:,:,1) = qXe;    PolaVecA(:,:,:,2) = qYe;
    PolaVecA(:,:,:,3) = qXe;    PolaVecA(:,:,:,4) = qYe;
elseif strcmp(p.polarization,'4EP')
    PolaVecA(:,:,:,1) = qXe;    PolaVecA(:,:,:,2) = qYe;
    PolaVecA(:,:,:,3) = qXre;   PolaVecA(:,:,:,4) = qYre;
elseif strcmp(p.polarization,'CHIDO')
    c = 1.2*pi; % from CHIDO paper
    PolaVecA(:,:,:,1) = cos(c*kr/p.kmax/2).*qL+1i*exp(-1i*phi).*sin(c*kr/p.kmax/2).*qR;
    PolaVecA(:,:,:,2) = cos(c*kr/p.kmax/2).*qR+1i*exp(1i*phi).*sin(c*kr/p.kmax/2).*qL;
    PolaVecA(:,:,:,3) = PolaVecA(:,:,:,1); PolaVecA(:,:,:,4) = PolaVecA(:,:,:,2);
elseif strcmp(p.polarization,'raMVR')
    qRa = cos(phi).*qX+sin(phi).*qY;
    qAz = -sin(phi).*qX+cos(phi).*qY;
    PolaVecA = zeros(p.Nk,p.Nk,3,p.Nc);
    PolaVecA(1:p.Nk/2,1:p.Nk/2,:,1) = qRa(1:p.Nk/2,1:p.Nk/2,:);
    PolaVecA(p.Nk/2+1:end,1:p.Nk/2,:,2) = qRa(p.Nk/2+1:end,1:p.Nk/2,:);
    for ii = 1:p.Nk/2
        PolaVecA(ii,ii:p.Nk-ii,:,3) = qAz(ii,ii:p.Nk-ii,:);
        PolaVecA(ii+1:end-ii+1,ii,:,4) = qAz(ii+1:end-ii+1,ii,:);
    end
end

% 4Pi channel setting
PolaVecB = flip(flip(PolaVecA,1),2);
if p.DualObj
    for ii=1:4
        PolaVecB(:,:,:,ii) = PolaVecB(:,:,:,ii).*exp(1i*cp(ii));
    end
else
    PolaVecB(:) = 0;
end

% aplanatic amplitude correction factor
Amplitude = pupilMask.*sqrt(cos_imm)./(n_med*cos_med);

% lateral translation
PhaseFactor = exp(1i*(kx*x0+ky*y0));


% pupil matraices of objective a b and derivative
PupilMatrixA = zeros(Nk,Nk,Nz,Nc,3);
PupilMatrixB = zeros(Nk,Nk,Nz,Nc,3);
dPupilMatrix = zeros(Nk,Nk,Nz,Nc,3,3);

for ii = 1:Nz, z = zl(ii);
    % propagation phase
    if strcmp(p.ztype,'medium')
        PhaseFactorA = PhaseFactor.*exp(1i*kz*z);
        PhaseFactorB = PhaseFactor.*exp(-1i*kz*z);
        % SAF
        if p.NA>p.n_med
            PhaseFactorA = PhaseFactorA.*exp(1i*(-p.fwd*p.kzimmnom+p.zstage*p.kzimm));
            PhaseFactorB = PhaseFactorB.*exp(1i*(-p.fwd*p.kzimmnom+p.zstage*p.kzimm)).*exp(2*1i*kz*p.depth);
        end
    elseif strcmp(p.ztype,'stage')
        PhaseFactorA = PhaseFactor.*exp(1i*kzimm*z);
        PhaseFactorB = PhaseFactor.*exp(-1i*kzimm*z);
    end

    for jj = 1:Nc
    for kk = 1:3
        PupilMatrixA(:,:,ii,jj,kk) = pupila.*Amplitude.*PhaseFactorA.*PolaVecA(:,:,kk,jj);
        PupilMatrixB(:,:,ii,jj,kk) = pupilb.*Amplitude.*PhaseFactorB.*PolaVecB(:,:,kk,jj);
    end
    end
end

% PSF engineering
% Vortex Phase Mask
if p.Vortex
    PupilMatrixA = PupilMatrixA.*exp(1i*phi);
    PupilMatrixB = PupilMatrixB.*exp(1i*phi);
end


% Pupil matrix (OTF) of the system, the composition of 2 objective's
PupilMatrix = PupilMatrixA+PupilMatrixB;

% normalized the total energy by in-focus free dipole
if p.NA>p.n_med
    PA0 = pupila.*Amplitude.*exp(1i*(-p.fwd*p.kzimmnom+p.zstage*p.kzimm)).*exp(1i*kz*p.depth).*PolaVecA;
    PB0 = pupilb.*Amplitude.*exp(1i*(-p.fwd*p.kzimmnom+p.zstage*p.kzimm)).*exp(1i*kz*p.depth).*PolaVecB;
else
    PA0 = pupila.*Amplitude.*PolaVecA;
    PB0 = pupilb.*Amplitude.*PolaVecB;
end
P0 = PA0+PB0;
Energy_norm = 1/3*sum(abs(P0).^2,"all");

% spatial derivatives of pupil matrix [dx,dy,dz]
dPupilMatrix(:,:,:,:,:,1) = 1i*kx.*PupilMatrix;       % dx
dPupilMatrix(:,:,:,:,:,2) = 1i*ky.*PupilMatrix;       % dy
if strcmp(p.ztype,'medium')  % dz
    dPupilMatrix(:,:,:,:,:,3) = 1i*kz.*(PupilMatrixA-PupilMatrixB);  
elseif strcmp(p.ztype,'stage')
    dPupilMatrix(:,:,:,:,:,3) = 1i*kzimm.*(PupilMatrixA-PupilMatrixB); 
end