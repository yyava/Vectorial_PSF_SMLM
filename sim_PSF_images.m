% This script display and output the images of Pupil functions,
% electric field components, and PSFs.

clc,clear
close all

p = set_parameters;

p.Nph = 4000;
p.Nbg = 12;

%% 
% p.detection =  "zslice"; 
% p.Nz = 1;
% p.zl = 0;   

p.polarization = 'CP';
p.DualObj = false;
p.Vortex = false;
% p.chi = 22.5;
p.Excitation = false;
p.dipoleType = 'diffusion';

p.azim = (45)/180*pi;
p.pola = 45/180*pi;
p.g2 = 0.0;

azimD = round(p.azim/pi*180);
polaD = round(p.pola/pi*180);
figTitle = strcat('\phi=',num2str(azimD),'° \theta=',num2str(polaD),'°');

p.z0 = 0e-9;

[PupilMatrix,dPupilMatrix,Energy_norm] = get_pupilMatrix(p);
[PSF,dPSF] = get_PSF(p,PupilMatrix,dPupilMatrix,Energy_norm);
M = p.Nx; w = p.czt_w; a = p.czt_a;
E = squeeze(czt2D(PupilMatrix,M,w,a)*(p.Dx*p.Dk/2/pi));
Energy = squeeze(sum(abs(PupilMatrix).^2,[1,2]));
[mu,dmu] = get_mu(p,PSF,dPSF); 

% 4-channel overview

plot_4Chan_2D(p,mu(:,:,ceil(p.Nz/2),:))
% colormap('hot')
sgtitle("4EP: "+figTitle+" z=0 nm NA=1.45",'FontSize',16,'FontWeight','bold')

% mu = 1e12*imnoise(mu*1e-12,'poisson');

[CRLBStore,~] = get_CRLB(p,mu,dmu);

%% PSF image: xy
ch = 1; % choose channel

T = PSF(:,:,ceil(p.Nz/2),ch)./max(PSF(:));
figure
imagesc([-1,1],[-1,1],T);
set(gca,'ydir','normal');
% colormap("hot");
clim([0,1])
axis equal
axis off
copygraphics(gca,'ContentType','vector')

%% PSF image: xz   % require zstack
T = squeeze(PSF(:,ceil(p.Nx/2),:,ch))./max(PSF(:));
figure
imagesc([-1,1],[-1,1],T');
hold on
set(gca,'ydir','normal');
colormap("hot");clim([0,1])
plot([-1,1],[0,0],'--g',LineWidth=1.5)
axis equal
axis off
copygraphics(gca,'ContentType','vector')

%% Pupil image: choose channel
PolarCh = 1;    % polarization component channel
EFCh = 1;       % electric field component channel
f_title = "P_{Xx}";         % Image name

Norm_med = p.n_med/p.NA;    % SAF range
T = PupilMatrix(:,:,1,PolarCh,EFCh)./max(abs(PupilMatrix(:)));
T(~p.pupilMask) = NaN;


%% Pupil image: complex amplitue
figure
% background
h1 = axes;
p1=imagesc([-1,1],[-1,1],zeros(size(T)));
set(p1,'alphadata',~isnan(T))
colormap(h1,"hot");clim([0,1])
axis equal
axis off
% title("$"+f_title+"$","Interpreter","latex")

% foreground
h2 = axes;
h=imagesc([-1,1],[-1,1],angle(T));
set(h,'alphadata',~isnan(T))
set(h,'alphadata',abs(T).^1.3)
set(gca,'ydir','normal');
colormap(h2,"hsv")
if Norm_med<=1
hold on
rectangle("Position",[-Norm_med,-Norm_med,2*Norm_med,2*Norm_med],"LineStyle","--","Curvature",[1,1],"LineWidth",1.5,"EdgeColor",'w')
end
clim([-pi,pi])
axis equal
axis off
% title("$"+f_title+"$","Interpreter","latex")
% fontsize(gcf,scale=2)

linkaxes([h1 h2])
copygraphics(gcf,'ContentType','vector')

%% Electric field image: choose channel
PolarCh = 1;    % polarization component channel
EFCh = 3;       % electric field component channel
f_title = "E_{Lx}"; % image name

T = E(:,:,PolarCh,EFCh)./max(abs(E(:)));
% phase correction of czt
T = T.*exp(1i*p.kmax*(p.xl)).*exp(1i*p.kmax*(p.xl'));

%% Electric field image: complex amplitude
figure
% background
h1 = axes;
p1=imagesc([-1,1],[-1,1],zeros(size(T)));
set(gca,'ydir','normal');
colormap(h1,"hot");clim([0,1])
axis equal
axis off
% title("$"+f_title+"$","Interpreter","latex")

% foreground
h2 = axes;
h=imagesc([-1,1],[-1,1],angle(T));
set(h,'alphadata',abs(T).^1.0)
colormap(h2,"hsv")
set(gca,'ydir','normal');
clim([-pi,pi])
axis equal
axis off
% title("$"+f_title+"$","Interpreter","latex")
% fontsize(gcf,scale=2)

linkaxes([h1 h2])
copygraphics(gcf,'ContentType','vector')

%% draw colorbar image for complex amplitude

CBa = ones(201,51);
CBp = ones(201,51);
CBa = CBa.*(0:50)/50;
CBp = CBp.*(0:200)'/200;

figure
% background
h1 = axes;
p1=imagesc([0,1],[0,2],zeros(size(CBa)));
set(gca,'ydir','normal');
colormap(h1,"hot");clim([0,1])
axis off
% foreground
h2 = axes;
h = imagesc([0,1],[0,2],CBp);
set(gca,'ydir','normal');
colormap(h2,"hsv")
set(h,'alphadata',CBa.^1.3)
axis off
linkaxes([h1 h2])

%% copy choosen image
copygraphics(gca,'ContentType','vector')


