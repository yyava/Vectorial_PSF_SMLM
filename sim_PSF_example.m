% This script gives an example how to use the functions in this file and
% how to simulate PSF in SMLM. Please run one section at a time.

clc,clear
% close all

p = set_parameters;

%% 1. Adjustable paramters
p.polarization = '4EP';
% polarization seeting: 'None','2LP','4LP','CP','2EP','4EP'
p.DualObj = false;      % 4Pi
p.Vortex = false;       % Vortex phase plate
p.Excitation = false;   % excitation modulation
p.z0 = 0e-9;            % axial position
p.Nph = 4000;           % signal photon counts
p.Nbg = 12;             % background photons per pixel across all channels

%% 2. simulation of 3D PSF
% free dipole
p.g2 = 0;

[PupilMatrix,dPupilMatrix,Energy_norm] = get_pupilMatrix(p);
[PSF_free,~] = get_PSF(p,PupilMatrix,dPupilMatrix,Energy_norm);

plot_4Chan_3D(p,PSF_free)
sgtitle('Free dipole PSF','FontSize',16,'FontWeight','bold')

%%
% fixed dipole 
p.g2 = 1.0;              % g2=1 for fieed dipole
p.azim = (45)/180*pi;    % polar angle
p.pola = 45.0/180*pi;    % azimuthal angle

[PupilMatrix,dPupilMatrix,Energy_norm] = get_pupilMatrix(p);    % pupil function
[PSF,dPSF] = get_PSF(p,PupilMatrix,dPupilMatrix,Energy_norm);   % PSF PDF

azimD = round(p.azim/pi*180);
polaD = round(p.pola/pi*180);
figTitle = strcat('\phi=',num2str(azimD),'° \theta=',num2str(polaD),'°');
plot_4Chan_3D(p,PSF)
sgtitle("Fixed dipole PSF: "+figTitle,'FontSize',16,'FontWeight','bold')

%% 3. Calculating CRLB with given parameters

[mu,dmu] = get_mu(p,PSF,dPSF);      % excepted photon counts
[CRLB,FisherM] = get_CRLB(p,mu,dmu);
CRLB(1:3,:) = CRLB(1:3,:)*1e9;      % meter to nm
CRLB(6:7,:) = CRLB(6:7,:)/pi*180;   % rad to degree
% Mean Angular Deviation
sigx = CRLB(6,:)*sin(p.pola);
sigy = CRLB(7,:);
a = (sigx.^2+sigy.^2)/2+abs((sigx.^2-sigy.^2)/2);
b = (sigx.^2+sigy.^2)/2-abs((sigx.^2-sigy.^2)/2);
[~,Eell] = ellipke(1-b./a);
CRLB(9,:,:) = sqrt(2/pi*a).*Eell;

Zl = p.zl*1e9;
figure
plot(Zl,CRLB(1,:),Zl,CRLB(2,:),Zl,CRLB(3,:),'LineWidth',1.5)
grid on
xlabel('z (nm)')
ylabel('Spatial CRLB (nm)')
legend('x','y','z'),legend(Box="off")
fontsize(gcf,scale=1.5)

figure
plot(Zl,CRLB(6,:),Zl,CRLB(7,:),Zl,CRLB(9,:),'LineWidth',1.5)
grid on
xlabel('z (nm)')
ylabel('Angular CRLB (°)')
legend('\phi','\theta','\delta'),legend(Box="off")
fontsize(gcf,scale=1.5)
