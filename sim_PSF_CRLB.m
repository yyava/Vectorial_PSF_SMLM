% This script plot the CRLB result by parameter sweeping

clc,clear
close all

p = set_parameters;
p.Nph = 4000;
p.Nbg = 12;

%% CRLB along z: sweeping
p.polarization = '4EP';
p.dipoleType = 'diffusion';
p.DualObj = false;
p.Vortex = false;
p.g2 = 0.75;
p.z0 = 0e-9;
p.Excitation = false;

% sample the polar and azimuthal angle

polal = acosd(0:0.05:0.95);
% aziml = 2.5:5:90-2.5;
% specific orientation
% polal = 45;
aziml = 45;

Ncfg = length(polal)*length(aziml);

Omega1 = zeros(1,Ncfg);
Omega2 = zeros(1,Ncfg);
for i = 1:length(aziml)
    for j = 1:length(polal)
        Omega1((i-1)*length(polal)+j) = aziml(i);
        Omega2((i-1)*length(polal)+j) = polal(j);
    end
end
CRLBStore = zeros(p.Np,p.Nz,Ncfg);

fprintf('\nStart fitting %i instances:\n',Ncfg); tic;

for jcfg = 1:Ncfg
    p.azim = Omega1(jcfg)/180*pi;
    p.pola = Omega2(jcfg)/180*pi;
    [PupilMatrix,dPupilMatrix] = get_pupilMatrix(p);
    [PSF,dPSF] = get_PSF(p,PupilMatrix,dPupilMatrix);
    [mu,dmu] = get_mu(p,PSF,dPSF); 
    [CRLBStore(:,:,jcfg),~] = get_CRLB(p,mu,dmu);
    disp(jcfg)
end
fprintf(['\nCRLB routine (cfg/second): ' num2str(toc,3) 's (' num2str(Ncfg/toc,5) ')\n'])

%% CRLB along z: data processing

CRLBStore(1:3,:,:) = CRLBStore(1:3,:,:)*1e9;
CRLBStore(6:7,:,:) = CRLBStore(6:7,:,:)/pi*180;
CRLBStore(6,:,:) = reshape(sind(Omega2),[1,1,Ncfg]).*CRLBStore(6,:,:);
% MAD
sigx = CRLBStore(6,:,:);
sigy = CRLBStore(7,:,:);
a = (sigx.^2+sigy.^2)/2+abs((sigx.^2-sigy.^2)/2);
b = (sigx.^2+sigy.^2)/2-abs((sigx.^2-sigy.^2)/2);
CRLBStore(9,:,:) = sqrt(pi)/8*a.*b.*(31*a.^2+31*b.^2+2*a.*b)./(a+b).^3.5;

CRLBz = zeros(4,p.Nz);
CRLBz(1,:) = mean((CRLBStore(1,:,:)+CRLBStore(2,:,:))/2,3);
CRLBz(2,:) = mean(CRLBStore(3,:,:),3);
CRLBz(3,:) = mean(CRLBStore(9,:,:),3);
CRLBz(4,:) = mean(CRLBStore(8,:,:),3);

% Store data for diff methods
% CRLBz_4EP = CRLBz;
% CRLBz_Vortex = CRLBz;
% CRLBz_CHIDO = CRLBz;
% CRLBz_raMVR = CRLBz;

%% CRLB along z: plot one method
Zl = p.zl*1e9;
% CRLB x z
figure("Position",[200,200,500,350])
plot(Zl,CRLBz(1,:),'b',Zl,CRLBz(2,:),'r',"LineWidth",1.5)
% ylim([0,4])
grid on
lgd = legend("$\sigma_x$ (nm)","$\sigma_z$ (nm)");
lgd.Interpreter = "latex";
lgd.Location = "northwest";
lgd.Box = "off";
xlabel('$z$ (nm)','Interpreter','latex')
ylabel('Localization precision','Interpreter','latex')
fontsize(gcf,scale=1.8)
% CRLB angle
CRLBz(5,:) = mean(CRLBStore(6,:,:),3);
CRLBz(6,:) = mean(CRLBStore(7,:,:),3);
figure("Position",[200,200,500,350])
plot(Zl,CRLBz(5,:),"Color",[.49,.73,.00],"LineWidth",1.5)
hold on
plot(Zl,CRLBz(6,:),"Color",[.00,.63,.95],"LineWidth",1.5)
plot(Zl,CRLBz(3,:),"Color",[.49,.18,.56],"LineWidth",1.5)
ylim([0,5])
grid on
lgd = legend("$\sin(\theta_d)\sigma_{\phi_d}\ (^{\circ})$", ...
    "$\sigma_{\theta_d}\ (^{\circ})$","$\bar{\delta}\ (^{\circ})$");
lgd.Interpreter = "latex";
lgd.Location = "northwest";
lgd.Box = "off";
xlabel('$z$ (nm)','Interpreter','latex')
ylabel('Angular precision','Interpreter','latex')
fontsize(gcf,scale=1.8)
% CRLB g2
figure("Position",[200,200,500,350])
plot(Zl,CRLBz(4,:),"Color",[1 0.5 0],"LineWidth",1.5)
ylim([0,0.1])
grid on
lgd = legend("$\sigma_{g_2}$");
lgd.Interpreter = "latex";
lgd.Location = "northwest";
lgd.Box = "off";
xlabel('$z$ (nm)','Interpreter','latex')
ylabel('$g_2$ precision','Interpreter','latex')
fontsize(gcf,scale=1.8)
%% CRLB along z: plot all methods
Zl = p.zl*1e9;
% CRLB x
figure("Position",[200,200,600,400])
% hold on
plot(Zl,CRLBz_Vortex(1,:),Zl,CRLBz_CHIDO(1,:), ...
    Zl,CRLBz_raMVR(1,:),Zl,CRLBz_4EP(1,:),"LineWidth",1.5)
ylim([0,12])
yticks(0:3:12)
grid on
lgd = legend("Vortex","CHIDO","raMVR","4EP");
lgd.Interpreter = "latex";
lgd.Location = "northwest";
lgd.Box = "off";
xlabel('$z$ (nm)','Interpreter','latex')
ylabel('$\sigma_x$ (nm)','Interpreter','latex')
fontsize(gcf,scale=1.8)
% CRLB z
figure("Position",[200,200,600,400])
plot(Zl,CRLBz_Vortex(2,:),Zl,CRLBz_CHIDO(2,:), ...
    Zl,CRLBz_raMVR(2,:),Zl,CRLBz_4EP(2,:),"LineWidth",1.5)
ylim([0,25])
yticks(0:5:25)
grid on
xlabel('$z$ (nm)','Interpreter','latex')
ylabel('$\sigma_z$ (nm)','Interpreter','latex')
fontsize(gcf,scale=1.8)
% CRLB delta
figure("Position",[200,200,600,400])
plot(Zl,CRLBz_Vortex(3,:),Zl,CRLBz_CHIDO(3,:), ...
    Zl,CRLBz_raMVR(3,:),Zl,CRLBz_4EP(3,:),"LineWidth",1.5)
ylim([0,8])
yticks(0:2:8)
grid on
xlabel('$z$ (nm)','Interpreter','latex')
ylabel('$\bar{\delta}\ (^{\circ})$','Interpreter','latex')
fontsize(gcf,scale=1.8)
% CRLB g2
figure("Position",[200,200,600,400])
plot(Zl,CRLBz_Vortex(4,:),Zl,CRLBz_CHIDO(4,:), ...
    Zl,CRLBz_raMVR(4,:),Zl,CRLBz_4EP(4,:),"LineWidth",1.5)

ylim([0,0.1])
yticks(0:0.02:0.1)
grid on
xlabel('$z$ (nm)','Interpreter','latex')
ylabel('$\sigma_{g_2}$','Interpreter','latex')
fontsize(gcf,scale=1.8)

copygraphics(gcf,'ContentType','vector')

%% CRLB w.r.t polar angle: sweeping
p.detection =  "zslice"; 
p.Nz = 1;
p.zl = 0;   

p.polarization = '4EP';
p.dipoleType = 'diffusion';
p.DualObj = false;
p.Vortex = false;
p.g2 = 0.75;
p.z0 = 0e-9;    % in focus
p.Excitation = false;

polal = acosd(0:0.01:1);
polal(end) = 1;     % avoid sigularity
aziml = 2.5:5:90-2.5;   % average
% aziml = 45;   % used one value when radial symmetry

CRLBStore = zeros(p.Np,length(polal),length(aziml));

for i=1:length(polal)
    for j=1:length(aziml)
        p.pola = polal(i)/180*pi;
        p.azim = aziml(j)/180*pi;
        [PupilMatrix,dPupilMatrix] = get_pupilMatrix(p);
        [PSF,dPSF] = get_PSF(p,PupilMatrix,dPupilMatrix);
        [mu,dmu] = get_mu(p,PSF,dPSF); 
        [CRLBStore(:,i,j),~] = get_CRLB(p,mu,dmu);
    end
    disp(i)
end

%% CRLB w.r.t polar angle: data processing

CRLBStore(1:3,:,:) = CRLBStore(1:3,:,:)*1e9;
CRLBStore(6:7,:,:) = CRLBStore(6:7,:,:)/pi*180;
CRLBStore(6,:,:) = sind(polal).*CRLBStore(6,:,:);
sigx = squeeze(CRLBStore(6,:,:));
sigy = squeeze(CRLBStore(7,:,:));
a = (sigx.^2+sigy.^2)/2+abs((sigx.^2-sigy.^2)/2);
b = (sigx.^2+sigy.^2)/2-abs((sigx.^2-sigy.^2)/2);
CRLBStore(9,:,:) = sqrt(pi)/8*a.*b.*(31*a.^2+31*b.^2+2*a.*b)./(a+b).^3.5;
CRLBStore = mean(CRLBStore,3);

% Store data for diff methods
% CRLBpola_4EP = CRLBStore;
% CRLBpola_Vortex = CRLBStore;
% CRLBpola_CHIDO = CRLBStore;
% CRLBpola_raMVR = CRLBStore;

%% CRLB w.r.t polar angle: plot one method
Xl = 1:length(polal);
% CRLB x z
figure("Position",[200,200,500,350])
plot(Xl,(CRLBStore(1,:)+CRLBStore(2,:))/2,'b', ...
    Xl,CRLBStore(3,:),'r',"LineWidth",1.5)
xlim([1,101])
xticks(cosd([90,75,60,45,30,15,1])*100+1)
xticklabels({'90','75','60','45','30','15  ','1'})
% ylim([0 4])
grid on
lgd = legend('$\sigma_x$ (nm)','$\sigma_z$ (nm)');
lgd.Interpreter = "latex";
lgd.Location = "north";
lgd.Box = "off";
xlabel('Polar angle $\theta_d\ (^{\circ})$','Interpreter','latex')
ylabel('Localization precision','Interpreter','latex')
fontsize(gcf,scale=1.8)

% CRLB angle
figure("Position",[200,200,500,350])
plot(Xl,CRLBStore(6,:),"Color",[.49,.73,.00],"LineWidth",1.5)
hold on
plot(Xl,CRLBStore(7,:),"Color",[.00,.63,.95],"LineWidth",1.5)
plot(Xl,CRLBStore(9,:),"Color",[.49,.18,.56],"LineWidth",1.5)
xlim([1,101])
xticks(cosd([90,75,60,45,30,15,1])*100+1)
xticklabels({'90','75','60','45','30','15  ','1'})
ylim([0,5])
grid on
lgd = legend("$\sin(\theta_d)\sigma_{\phi_d}\ (^{\circ})$", ...
    "$\sigma_{\theta_d}\ (^{\circ})$","$\bar{\delta}\ (^{\circ})$");
lgd.Interpreter = "latex";
lgd.Location = "northwest";
lgd.Box = "off";
xlabel('Polar angle $\theta_d\ (^{\circ})$','Interpreter','latex')
ylabel('Angular precision','Interpreter','latex')
fontsize(gcf,scale=1.8)
% CRLB g2
figure("Position",[200,200,500,350])
plot(Xl,CRLBStore(8,:),"Color",[1 0.5 0],"LineWidth",1.5)
xlim([1,101])
xticks(cosd([90,75,60,45,30,15,1])*100+1)
xticklabels({'90','75','60','45','30','15  ','1'})
ylim([0,0.1])
grid on
lgd = legend("$\sigma_{g_2}$");
lgd.Interpreter = "latex";
lgd.Location = "northwest";
lgd.Box = "off";
xlabel('Polar angle $\theta_d\ (^{\circ})$','Interpreter','latex')
ylabel('$g_2$ precision','Interpreter','latex')
fontsize(gcf,scale=1.8)

%% CRLB w.r.t polar angle: plot all methods
Xl = 1:length(polal);
% CRLB x
figure("Position",[200,200,600,400])
plot(Xl,CRLBpola_Vortex(1,:),Xl,CRLBpola_CHIDO(1,:), ...
    Xl,CRLBpola_raMVR(1,:),Xl,CRLBpola_4EP(1,:),"LineWidth",1.5)
xlim([1,101])
xticks(cosd([90,75,60,45,30,15,1])*100+1)
xticklabels({'90','75','60','45','30','15  ','1'})
ylim([0,12])
yticks(0:3:12)
grid on
lgd = legend("Vortex","CHIDO","raMVR","4EP");
lgd.Interpreter = "latex";
lgd.Location = "northwest";
lgd.Box = "off";
xlabel('Polar angle $\theta_d\ (^{\circ})$','Interpreter','latex')
ylabel('$\sigma_x$ (nm)','Interpreter','latex')
fontsize(gcf,scale=1.8)
% CRLB delta
figure("Position",[200,200,600,400])
plot(Xl,CRLBpola_Vortex(6,:),"Color",[.00 .45 .74],"LineWidth",1.5)
hold on 
plot(Xl,CRLBpola_CHIDO(6,:),"Color",[.85 .33 .10],"LineWidth",1.5)
plot(Xl,CRLBpola_raMVR(6,:),"Color",[.93 .69 .13],"LineWidth",1.5)
plot(Xl,CRLBpola_4EP(6,:),"Color",[.49 .18 .56],"LineWidth",1.5)
xlim([1,101])
xticks(cosd([90,75,60,45,30,15,1])*100+1)
xticklabels({'90','75','60','45','30','15  ','1'})
ylim([0,10])
grid on
xlabel('Polar angle $\theta_d\ (^{\circ})$','Interpreter','latex')
ylabel('$\sin(\theta_d)\sigma_{\phi_d}\ (^{\circ})$','Interpreter','latex')
fontsize(gcf,scale=1.8)

plot(Xl,CRLBpola_Vortex(7,:),"LineStyle",'--',"Color",[.00 .45 .74],"LineWidth",1.5)
plot(Xl,CRLBpola_CHIDO(7,:),"LineStyle",'--',"Color",[.85 .33 .10],"LineWidth",1.5)
plot(Xl,CRLBpola_raMVR(7,:),"LineStyle",'--',"Color",[.93 .69 .13],"LineWidth",1.5)
plot(Xl,CRLBpola_4EP(7,:),"LineStyle",'--',"Color",[.49 .18 .56],"LineWidth",1.5)
plot(Xl,CRLBpola_Vortex(9,:),"Color",[.00 .45 .74],"LineWidth",1.5)
plot(Xl,CRLBpola_CHIDO(9,:),"Color",[.85 .33 .10],"LineWidth",1.5)
plot(Xl,CRLBpola_raMVR(9,:),"Color",[.93 .69 .13],"LineWidth",1.5)
plot(Xl,CRLBpola_4EP(9,:),"Color",[.49 .18 .56],"LineWidth",1.5)

h1 = plot(-1,0,'k:');
h2 = plot(-1,0,'k--');
h3 = plot(-1,0,'k');
xlim([1,101])
xticks(cosd([90,75,60,45,30,15,1])*100+1)
xticklabels({'90','75','60','45','30','15  ','1'})
ylim([0,20])
grid on
lgd = legend([h1,h2,h3],'$\sin(\theta_d)\sigma_{\phi_d}\ (^{\circ})$', ...
    '$\sigma_{\theta_d}\ (^{\circ})$','$\bar{\delta}\ (^{\circ})$');
lgd.Interpreter = "latex";
lgd.Location = "north";
lgd.Box = "off";
xlabel('Polar angle $\theta_d\ (^{\circ})$','Interpreter','latex')
ylabel('$\sigma_x$ (nm)','Interpreter','latex')
fontsize(gcf,scale=1.8)

%% CRLB w.r.t azimuthal angle: sweeping
p.detection =  "zslice"; 
p.Nz = 1;
p.zl = 0;   

p.polarization = '4EP';
p.dipoleType = 'diffusion';
p.DualObj = false;
p.Vortex = false;
p.g2 = 0.75;
p.z0 = 0e-9;    % in focus
p.Excitation = false;

aziml = 0:1:180;
polal = acosd(0:0.05:0.95);
CRLBStore = zeros(p.Np,length(polal),length(aziml));

for i=1:length(polal)
    for j=1:length(aziml)
        p.pola = polal(i)/180*pi;
        p.azim = aziml(j)/180*pi;
        [PupilMatrix,dPupilMatrix] = get_pupilMatrix(p);
        [PSF,dPSF] = get_PSF(p,PupilMatrix,dPupilMatrix);
        [mu,dmu] = get_mu(p,PSF,dPSF); 
        [CRLBStore(:,i,j),~] = get_CRLB(p,mu,dmu);
    end
    disp(i)
end

%% CRLB w.r.t azimuthal angle: data processing
CRLBStore(1:3,:,:) = CRLBStore(1:3,:,:)*1e9;
CRLBStore(6:7,:,:) = CRLBStore(6:7,:,:)/pi*180;
CRLBStore(6,:,:) = sind(polal).*CRLBStore(6,:,:);
sigx = CRLBStore(6,:,:);
sigy = CRLBStore(7,:,:);
a = (sigx.^2+sigy.^2)/2+abs((sigx.^2-sigy.^2)/2);
b = (sigx.^2+sigy.^2)/2-abs((sigx.^2-sigy.^2)/2);
CRLBStore(9,:,:) = sqrt(pi)/8*a.*b.*(31*a.^2+31*b.^2+2*a.*b)./(a+b).^3.5;
CRLBStore = squeeze(mean(CRLBStore,2));

%% CRLB w.r.t azimuthal angle: plot one method
% CRLB x
figure("Position",[200,200,500,350])
plot(aziml,(CRLBStore(1,:)+CRLBStore(2,:))/2,'b',"LineWidth",1.5)
xlim([1,180])
xticks(0:45:180)
ylim([0 4])
grid on
lgd = legend('$\sigma_x$ (nm)');
lgd.Interpreter = "latex";
lgd.Location = "north";
lgd.Box = "off";
xlabel('Azimuthal angle $\phi_d\ (^{\circ})$','Interpreter','latex')
ylabel('Lateral precision','Interpreter','latex')
fontsize(gcf,scale=1.8)
% CRLB z
figure("Position",[200,200,500,350])
plot(aziml,CRLBStore(3,:),'r',"LineWidth",1.5)
xlim([1,180])
xticks(0:45:180)
ylim([0 20])
grid on
lgd = legend('$\sigma_z$ (nm)');
lgd.Interpreter = "latex";
lgd.Location = "north";
lgd.Box = "off";
xlabel('Azimuthal angle $\phi_d\ (^{\circ})$','Interpreter','latex')
ylabel('Lateral precision','Interpreter','latex')
fontsize(gcf,scale=1.8)
% CRLB angle
figure("Position",[200,200,500,350])
plot(aziml,CRLBStore(6,:),"Color",[.49,.73,.00],"LineWidth",1.5)
hold on
plot(aziml,CRLBStore(7,:),"Color",[.00,.63,.95],"LineWidth",1.5)
plot(aziml,CRLBStore(9,:),"Color",[.49,.18,.56],"LineWidth",1.5)
xlim([1,180])
xticks(0:45:180)
ylim([0 5])
grid on
lgd = legend('$\sin(\theta_d)\sigma_{\phi_d}\ (^{\circ})$', ...
    '$\sigma_{\theta_d}\ (^{\circ})$','$\bar{\delta}\ (^{\circ})$');
lgd.Interpreter = "latex";
lgd.Location = "north";
lgd.Box = "off";
xlabel('Azimuthal angle $\phi_d\ (^{\circ})$','Interpreter','latex')
ylabel('Angular precision','Interpreter','latex')
fontsize(gcf,scale=1.8)
%CRLB g2
figure("Position",[200,200,500,350])
plot(aziml,CRLBStore(8,:),"Color",[1 0.5 0],"LineWidth",1.5)
xlim([1,180])
xticks(0:45:180)
ylim([0 0.06])
grid on
lgd = legend('$\sigma_{g_2}$');
lgd.Interpreter = "latex";
lgd.Location = "north";
lgd.Box = "off";
xlabel('Azimuthal angle $\phi_d\ (^{\circ})$','Interpreter','latex')
ylabel('$g_2$ precision','Interpreter','latex')
fontsize(gcf,scale=1.8)

%% CRLB w.r.t g2: sweeping
p.detection =  "zslice"; 
p.Nz = 1;
p.zl = 0;   

p.polarization = '4EP';
p.dipoleType = 'diffusion';
p.DualObj = false;
p.Vortex = false;
p.g2 = 0.75;
p.z0 = 0e-9;    % in focus
p.Excitation = false;

polal = acosd(0:0.05:0.95);
aziml = 2.5:5:90-2.5;
g2l = 0:0.025:1;

Ncfg = length(polal)*length(aziml);
Omega1 = zeros(1,Ncfg);
Omega2 = zeros(1,Ncfg);
for i = 1:length(aziml)
    for j = 1:length(polal)
        Omega1((i-1)*length(polal)+j) = aziml(i);
        Omega2((i-1)*length(polal)+j) = polal(j);
    end
end
CRLBStore = zeros(p.Np,length(g2l),Ncfg);
for i = 1:length(g2l)
    p.g2 = g2l(i);
    for jcfg = 1:Ncfg
    p.azim = Omega1(jcfg)/180*pi;
    p.pola = Omega2(jcfg)/180*pi;
    [PupilMatrix,dPupilMatrix] = get_pupilMatrix(p);
    [PSF,dPSF] = get_PSF(p,PupilMatrix,dPupilMatrix);
    [mu,dmu] = get_mu(p,PSF,dPSF); 
    [CRLBStore(:,i,jcfg),~] = get_CRLB(p,mu,dmu);
    end
    disp(i)
end

%% CRLB w.r.t g2: data processing
CRLBStore(1:3,:,:) = CRLBStore(1:3,:,:)*1e9;
CRLBStore(6:7,:,:) = CRLBStore(6:7,:,:)/pi*180;
CRLBStore(6,:,:) = reshape(sind(Omega2),[1,1,Ncfg]).*CRLBStore(6,:,:);
sigx = CRLBStore(6,:,:);
sigy = CRLBStore(7,:,:);
a = (sigx.^2+sigy.^2)/2+abs((sigx.^2-sigy.^2)/2);
b = (sigx.^2+sigy.^2)/2-abs((sigx.^2-sigy.^2)/2);
CRLBStore(9,:,:) = sqrt(pi)/8*a.*b.*(31*a.^2+31*b.^2+2*a.*b)./(a+b).^3.5;
CRLBStore = squeeze(mean(CRLBStore,3));

%% CRLB w.r.t g2: plot
% CRLBx
figure("Position",[200,200,500,350])
plot(g2l,(CRLBStore(1,:)+CRLBStore(2,:))/2,'b',"LineWidth",1.5)
% xlim([1,180])
% xticks(0:45:180)
ylim([0 4])
grid on
lgd = legend('$\sigma_x$ (nm)');
lgd.Interpreter = "latex";
lgd.Location = "north";
lgd.Box = "off";
xlabel('Orientation constraint $g_2$','Interpreter','latex')
ylabel('Lateral precision','Interpreter','latex')
fontsize(gcf,scale=1.8)
% CRLBz
figure("Position",[200,200,500,350])
plot(g2l,CRLBStore(3,:),'r',"LineWidth",1.5)
xlim([0,1])
% xticks(0:45:180)
ylim([0 20])
grid on
lgd = legend('$\sigma_z$ (nm)');
lgd.Interpreter = "latex";
lgd.Location = "north";
lgd.Box = "off";
xlabel('Orientation constraint $g_2$','Interpreter','latex')
ylabel('Lateral precision','Interpreter','latex')
fontsize(gcf,scale=1.8)
% CRLBd
figure("Position",[200,200,500,350])
plot(g2l,CRLBStore(6,:),'c',g2l,CRLBStore(7,:),'g', ...
    g2l,CRLBStore(9,:),'m',"LineWidth",1.5)
xlim([0,1])
% xticks(0:45:180)
ylim([0 4])
grid on
lgd = legend('$\sin(\theta_d)\sigma_{\phi_d}\ (^{\circ})$', ...
    '$\sigma_{\theta_d}\ (^{\circ})$','$\bar{\delta}\ (^{\circ})$');
lgd.Interpreter = "latex";
lgd.Location = "north";
lgd.Box = "off";
xlabel('Orientation constraint $g_2$','Interpreter','latex')
ylabel('Angular precision','Interpreter','latex')
fontsize(gcf,scale=1.8)
%CRLBg2
figure("Position",[200,200,500,350])
plot(g2l,CRLBStore(8,:),"Color",[1 0.5 0],"LineWidth",1.5)
xlim([0,1])
% xticks(0:45:180)
ylim([0 0.06])
grid on
lgd = legend('$\sigma_{g_2}$');
lgd.Interpreter = "latex";
lgd.Location = "north";
lgd.Box = "off";
xlabel('Orientation constraint $g_2$','Interpreter','latex')
ylabel('$g_2$ precision','Interpreter','latex')
fontsize(gcf,scale=1.8)







