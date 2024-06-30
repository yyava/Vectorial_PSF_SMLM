function p = set_parameters
% This function sets all parameters for vectorial PSF calculations
% @Menglin Wu

%% optical parameters
p.NA = 1.45;        % numerical aperture
p.n_med = 1.333;    % refractive indices of medium
p.n_cov = 1.523;    % refractive indices of coverslip
p.n_imm = 1.518;    % refractive indices of immersion environment
% oil: 1.518; air: 1.0
% nominal refractive index of immersion fluid matching objective lens design
p.n_immnom = p.n_cov;   
p.fwd = 120e-6;    % nominal working distance of objective lens
p.zstage = 118.935e-6;  % (NA=1.45, λ=520, f=200, nomf=700)
p.depth = 200e-9;        % distance of image plane from coverslip
p.zrange = [-500e-9,500e-9];    % axial range in fitting [min,max]
% flag for axial range by z-position in medium or by z-stage position
p.ztype = 'medium'; %'stage';

% p.lambda = 597.5e-9;   % wavelength
p.lambda = 520e-9;   % wavelength
p.k0 = 2*pi/p.lambda;   % wavenumber in air

% emitter position (shift inside medium)
p.x0 = 0;
p.y0 = 0;
p.z0 = 0;

%% simulation sampling
p.magObj = 100;      % magnification of objective
p.cameraSize = 6.5e-6;  % pixel size of camera is 6.5 um

% Variable naming rules:
% D:pixel pitch; N: sampling number; L: window width; l: 1D array
% Example: Lx=Dx*Nx; xl=0:Dx:Lx-Dx;

% Image plane: x, y direction (same pixel size in x and y)
p.Dx = p.cameraSize/p.magObj;   % image pixel size
p.Nx = 19;
% smaller pixel size for plotting PSF images
% p.Dx = 5e-9;
% p.Nx = 201;

p.Lx = p.Nx*p.Dx;   % length of ROI
p.xl = -(p.Lx-p.Dx)/2:p.Dx:(p.Lx-p.Dx)/2;   % coordinate array

% Fourier plane
p.kmax = p.NA*p.k0;     % maximum angular frquency
p.Nk = 32;      % sampling number in pupil plane
% smaller pixel size for plotting pupil images
% p.Nk = 256;   % 400 for plotting pupil function
p.Dk = p.kmax*2/p.Nk;
p.kl = -(p.Nk-1)/2*p.Dk:p.Dk:(p.Nk-1)/2*p.Dk;    % size:(1,Nk)
[p.kx, p.ky] = meshgrid(p.kl,p.kl); % size:(Nk,Nk)
p.kr = sqrt(p.kx.^2+p.ky.^2);       % size:(Nk,Nk)
p.kz = sqrt(p.n_med^2*p.k0^2-p.kr.^2);      % wave vector z; size:(Nk,Nk)
p.kzimm = sqrt(p.n_imm^2*p.k0^2-p.kr.^2);   % wave vector in immersion fluid
p.kzimmnom = sqrt(p.n_immnom^2*p.k0^2-p.kr.^2); 
p.phi = atan2(p.ky,p.kx);           % azimuth; size:(Nk,Nk)

% z direction
% zslice only contains emitter plane; zstack contains all z-position
p.detection =  "zstack"; %"zslice"; 

p.Lz = 400e-9;
p.Dz = 5e-9;
p.Nz = round(p.Lz/p.Dz+1);
p.zl = -p.Lz/2:p.Dz:p.Lz/2;
% % p.zl = 0:p.Dz:p.Lz; % half
% p.zl = -200e-9:p.Dz:400e-9;

if p.detection == "zslice"
    p.Nz = 1;
    p.zl = 0;
end

%% default dipole parameters setting
% when dipole type is fixed or free, g2 will be ignored
d = 'diffusion'; f = 'fixed'; r = 'free';
p.dipoleType = d;
p.pola = 60.0/180*pi;   % polar angle
p.azim = 45.0/180*pi;   % azimuthal angle
p.g2 = 0.75;               % g2=1:fixed; g2=0:free

%% channel setting
% polarization seeting: 'None','2LP','4LP','CP','2EP','4EP'
p.polarization = '2LP';
p.chi = 22.5;   % rotation angle of QWP for EP

p.Vortex = false;   % Vortex Phase

% 4Pi setting
p.DualObj = false;  % 4Pi 

p.ta = 1.;      % transmission coefficient
p.tb = 1.;
p.phi0 = 0;     % phi_a-phi_b, phase diff of 2 arms
p.phixy = pi/2; % phi_x-phi_y, phase diff between x and y polarization
p.Nc = 4;       % number of channels
p.chPhase = [p.phi0,p.phi0+pi,...
    p.phi0+p.phixy,p.phi0+p.phixy+pi];  % channel phase

%% excitation modulation setting
p.Excitation = false;
p.M = 3;                % number of excitation states
p.Xi = [0,pi/3,pi*2/3]; % linear polarization angle
% p.M = 4;
% p.Xi = [0,pi/4,pi/2,pi*3/4];
%% photon count
p.Nph = 3000;   % signal photon count
p.Nbg = 4;     % background photon count per pixel

%% aberration setting
% Phase of pupil plane
p.Wa = zeros(p.Nk);
p.Wb = zeros(p.Nk);

% Zernike orders [n1,m1,A1;n2,m2,A2;...] 
% with n1,n2,... the radial orders, 
% m1,m2,... the azimuthal orders, 
% and A1,A2,... the Zernike coefficients in lambda rms, 
% A = 0.072 means diffraction limit
p.aberrations = [ ...
    2, -2,  0;  % Oblique astigmatism
    2,  2,  0;  % Vertical astigmatism
    3, -1,  0;  % Vertical coma
    3,  1,  0;  % Horizontal coma
    3, -3,  0;  % Vertical trefoil
    3,  3,  0;  % Oblique trefoil
    4,  0,  0;  % Primary spherical
    4, -2,  0;
    4,  2,  0;
    5, -1,  0;
    5,  1,  0;
    6,  0,  0];
p.Wa = get_ZernikeFunc(p.Nk,p.aberrations);
p.Wb = get_ZernikeFunc(p.Nk,p.aberrations);

%% pupil setting
p.pupilMask = (p.kr<p.kmax);   % pupil window shape
p.pupila = p.ta*p.pupilMask.*exp(1i*p.Wa); % coherent OTF of obj.a
p.pupilb = p.tb*p.pupilMask.*exp(1i*p.Wb); % coherent OTF of obj.b

%% fiting model
p.fitModel = 'xyz-azim-pola-diffusion'; % difussion dipole
% p.fitModel = 'xyz-azim-pola'; % fixed dipole
% p.fitModel = 'xyz';   % free dipole
% p.fitModel = 'xy-azim-pola-diffusion';
if      contains(p.fitModel,'diffusion')
    % Theta = [dx,dy,dz,Nph,Nbg,azim,pola,g2]
    p.Np = 8;
elseif  contains(p.fitModel,'azim-pola')
    % Theta = [dx,dy,dz,Nph,Nbg,azim,pola]
    p.Np = 7;
else
    % Theta = [dx,dy,dz,Nph,Nbg]
    p.Np = 5;
end

if contains(p.fitModel,'diffusion')
    if ~strcmp(p.dipoleType,'diffusion')
        error("Fitting diffusion coefficient is not supported for fixed or free dipole...\n" + ...
            "Change the dipole type")
    end
end

%% calculate additional vectors for chirp z transform
% parameters of MATLAB czt function
p.czt_w = exp(-1i*p.Dk*p.Dx);   % Ratio between spiral contour points
p.czt_a = exp(-1i*(p.Lx-p.Dx)/2*p.Dk);  % Spiral contour initial point

%% fitting parameters in localization
p.NiterMax = 30;    % max iteration number
p.tollim = 1e-6;    % tolerent limit
p.varfit = 0;       % Gaussian readout noise
p.Ncfg = 1;         % number of configure
p.flg_parallel = false; % parallel computing flag

