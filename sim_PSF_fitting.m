% Simulate and localize dipole emitter

clc,clear
close all

p = set_parameters;

p.Ncfg = 10000;

p.polarization = 'CP';
p.dipoleType = 'fixed';
p.DualObj = false;
p.Excitation = false;

p.detection =  "zslice"; 
p.Nz = 1;
p.zl = 0; 

%% generate ground truth PSFs

allmu = zeros(p.Nx,p.Nx,p.Nz,p.Nc,p.Ncfg);
alldmu= zeros(p.Nx,p.Nx,p.Nz,p.Nc,p.Np,p.Ncfg);
object = zeros(p.Np,p.Ncfg);
if p.Excitation
    allmu = zeros(p.Nx,p.Nx,p.Nz,p.Nc*p.M,p.Ncfg);
    alldmu= zeros(p.Nx,p.Nx,p.Nz,p.Nc*p.M,p.Np,p.Ncfg);
end

for ii = 1:p.Ncfg
    % true parameters
    % dx = (1-2*rand)*p.Dx;
    % dy = (1-2*rand)*p.Dx;
    % dz = (1-2*rand)*400e-9;
    % dazim = rand*2*pi;  % 0~2pi
    % dpola = acos(rand); % 0~pi/2
    % dg2 = 0.8+(1-2*rand)*0.2;
    dx = 0;
    dy = 0;
    dz = 0e-9;
    Nph = 4000;
    Nbg = 12;
    dazim = (45)/180*pi;  % 0~2pi
    dpola = 45.0/180*pi; % 0~pi/2
    dg2 = 0.8;

    % generate object and PSFs based on fitModel
    % object(:,ii) = [dx,dy,dz,Nph,Nbg];
    % object(:,ii) = [dx,dy,dz,Nph,Nbg,dazim,dpola];
    object(:,ii) = [dx,dy,dz,Nph,Nbg,dazim,dpola,dg2];
    [allmu(:,:,:,:,ii),alldmu(:,:,:,:,:,ii)] = get_PoissonRate(p,object(:,ii));
end

% add noise
allspots = 1e12*imnoise(allmu*1e-12,'poisson');

%% %% calculate CRLB for each set of parameter
CRLBstore = zeros(p.Np,p.Ncfg);
for i = 1:p.Ncfg
    [CRLBstore(:,i),~] = get_CRLB(p,allmu(:,:,:,:,i),alldmu(:,:,:,:,:,i));
end

CRLBstore(6:7,:) = CRLBstore(6:7,:)/pi*180;
CRLBstore(1:3,:) = CRLBstore(1:3,:)*1e9;  
sigx = sin(dpola)*CRLBstore(6,:);
sigy = CRLBstore(7,:);
a = (sigx.^2+sigy.^2)/2+abs((sigx.^2-sigy.^2)/2);
b = (sigx.^2+sigy.^2)/2-abs((sigx.^2-sigy.^2)/2);
% CRLB(9,:) = sqrt(pi)/8*a.*b.*(31*a.^2+31*b.^2+2*a.*b)./(a+b).^3.5;
[~,Eell] = ellipke(1-b./a);
CRLBstore(9,:) = sqrt(2/pi*a).*Eell;

CRLBmean = mean(CRLBstore,2);
%% MLE fitting routine
Theta0 = initValues(allspots,p);
% Theta0 = object;
p.flg_parallel = false; % parallel computing

[ThetaStore,muStore,dmuStore,meritStore,numiters,ThetaAllStore] = localization(allspots,Theta0,p);

%% Fitting PSF image overview
cfg = 1;
ExCh = 1;   % excitation modulation channel
plot_4Chan_2D(p,allmu(:,:,1,4*(ExCh-1)+1:4*ExCh,cfg))
sgtitle("True image")
plot_4Chan_2D(p,allspots(:,:,1,4*(ExCh-1)+1:4*ExCh,cfg))
sgtitle("Noisy image")
[allmufit,~] = get_PoissonRate(p,ThetaStore(:,cfg));
plot_4Chan_2D(p,allmufit(:,:,1,4*(ExCh-1)+1:4*ExCh))
sgtitle("Fitting image")

%% Statistic
outliners = (numiters==p.NiterMax+1);
%%%% remove unconvergent result
% ThetaStore = ThetaStore(:,~outliners);
% object = object(:,~outliners);

DTheta = ThetaStore-object; % fitting deviation
DTheta(1:3,:) = DTheta(1:3,:)*1e9;
DTheta(6:7,:) = DTheta(6:7,:)/pi*180;
DTheta(6,:) = mod(DTheta(6,:)+180,360)-180;
%%%% angle between estimated and true orientation
    v1x = sin(ThetaStore(7,:)).*cos(ThetaStore(6,:));
    v2x = sin(object(7,:)).*cos(object(6,:));
    v1y = sin(ThetaStore(7,:)).*sin(ThetaStore(6,:));
    v2y = sin(object(7,:)).*sin(object(6,:));
    v1z = cos(ThetaStore(7,:));
    v2z = cos(object(7,:));
DTheta(9,:) = acos(abs(v1x.*v2x+v1y.*v2y+v1z.*v2z));
DTheta(9,:) = DTheta(9,:)/pi*180;

DThetaMean = mean(DTheta,2);
DThetaStd = std(DTheta,0,2);

%% Plot results
%%%% 3D localization result
figure("Position",[100,100,600,600])
S = 3;
C = DTheta(3,:);
scatter3(DTheta(1,:),DTheta(2,:),DTheta(3,:),S,C,"filled")
hold on
[X,Y,Z] = ellipsoid(DThetaMean(1),DThetaMean(2),DThetaMean(3),DThetaStd(1),DThetaStd(2),DThetaStd(3));
surf(X,Y,Z,"FaceAlpha",0.3,"EdgeColor","none","FaceColor","red");
% [X,Y,Z] = ellipsoid(0,0,0,CRLB(1),CRLB(2),CRLB(3));
% surf(X,Y,Z,"FaceAlpha",0.3,"EdgeColor","none","FaceColor","blue");
xlim([-10,10])
ylim([-10,10])
zlim([-10,10])
xlabel("$\Delta x$ (nm)",'Interpreter','latex')
ylabel("$\Delta y$ (nm)",'Interpreter','latex')
zlabel("$\Delta z$ (nm)",'Interpreter','latex')
fontsize(gcf,scale=1.8)
% title("z0=50 nm, \phi=90°, \theta=80°")

copygraphics(gcf,'ContentType','vector')

%%% 2D localization result
figure("Position",[100,100,500,500])
S = 3;
scatter(DTheta(1,:),DTheta(2,:),S,'k',"filled")
hold on
rectangle('Position', [DThetaMean(1)-DThetaStd(1),DThetaMean(2)-DThetaStd(2),2*DThetaStd(1),2*DThetaStd(2)], 'Curvature', [1 1],'EdgeColor', 'r');
axis equal
grid on
box on
xlim([-10,10])
ylim([-10,10])
xlabel("$\Delta x$ (nm)",'Interpreter','latex')
ylabel("$\Delta y$ (nm)",'Interpreter','latex')
fontsize(gcf,scale=1.8)

%% Angular deviation
% figure("Position",[200,200,350,350])
% axes('Position',[0.25, 0.2, 0.7, 0.75])
% figure("Position",[200,200,350,350])
% axes('Position',[0.25, 0.2, 0.7, 0.75])
figure()
histogram(DTheta(9,:),30,Normalization="pdf",BinLimits=[0,6])
% histogram(DTheta(9,:),50,BinLimits=[0,10])
hold on 
t1=plot([DThetaMean(9),DThetaMean(9)],[0,0.6],'b',"LineWidth",2,"DisplayName","Simulation");
plot([CRLBmean(9),CRLBmean(9)],[0,0.6],'r--',"LineWidth",2)
xlim([0 6])
xticks(0:1:6)
grid on
% ylim([0,0.35])
t = 0:0.01:6;
a=a(1);b=b(1);
Gd = t/sqrt(a*b).*exp(-t.^2/4*(1/a+1/b)).*besseli(0,t.^2/4*(1/b-1/a));
% Gd = t/a/b.*exp(-t.^2/4*(1/a^2+1/b^2)).*(1+t.^4/64*(1/a^2-1/b^2)^2);
t2=plot(t,Gd,'r',LineWidth=2,DisplayName="Theory");
legend([t1,t2])
xlabel("$\delta (^\circ)$",'Interpreter','latex')
ylabel('Probability density','Interpreter','latex')
fontsize(gcf,scale=1.8)

%% g2
t = -0.2:0.001:0.2;
ft = normpdf(t,DThetaMean(8),DThetaStd(8));

figure("Position",[200,200,350,350])
axes('Position',[0.25, 0.2, 0.7, 0.75])
histogram(DTheta(8,:),41,BinLimits=[-0.2,0.2],Normalization="pdf",FaceColor=[.85,.33,.10])
hold on
plot(t,ft,'k',"LineWidth",1.5)
plot([DThetaMean(8),DThetaMean(8)],[0,15],'r',"LineWidth",1.5)
xlim([-0.2,0.2])
ylim([0,10])
grid on
xlabel("$\Delta g_2$",'Interpreter','latex')
ylabel('Probability density','Interpreter','latex')
fontsize(gcf,scale=1.8)

%% x
t = -10:0.1:10;
ft = normpdf(t,DThetaMean(1),DThetaStd(1));
figure("Position",[200,200,350,350])
axes('Position',[0.25, 0.2, 0.7, 0.75])
histogram(DTheta(1,:),41,Normalization="pdf",BinLimits=[-10,10],FaceColor=[.49,.73,.00])
hold on
plot(t,ft,'k',"LineWidth",1.5)
plot([DThetaMean(1),DThetaMean(1)],[0,1],'r',"LineWidth",1.5)
xlim([-10,10])
ylim([0,0.3])
grid on
xlabel("$\Delta x$ (nm)",'Interpreter','latex')
ylabel('Probability density','Interpreter','latex')
fontsize(gcf,scale=1.8)
%% y
ft = normpdf(t,DThetaMean(2),DThetaStd(2));
figure("Position",[200,200,350,350])
axes('Position',[0.25, 0.2, 0.7, 0.75])
histogram(DTheta(2,:),41,Normalization="pdf",BinLimits=[-10,10],FaceColor=[.00,.63,.95])
hold on
plot(t,ft,'k',"LineWidth",1.5)
plot([DThetaMean(2),DThetaMean(2)],[0,1],'r',"LineWidth",1.5)
xlim([-10,10])
ylim([0,0.3])
grid on
xlabel("$\Delta y$ (nm)",'Interpreter','latex')
ylabel('Probability density','Interpreter','latex')
fontsize(gcf,scale=1.8)
%% z
len = 4;
t = -len:0.1:len;
ft = normpdf(t,DThetaMean(3),DThetaStd(3));
figure("Position",[200,200,350,350])
axes('Position',[0.25, 0.2, 0.7, 0.75])
histogram(DTheta(3,:),41,Normalization="pdf",BinLimits=[-len,len],FaceColor=[.49,.18,.56])
hold on
plot(t,ft,'k',"LineWidth",1.5)
plot([DThetaMean(3),DThetaMean(3)],[0,1],'r',"LineWidth",1.5)
xlim([-len,len])
ylim([0,0.8])
grid on
xlabel("$\Delta z$ (nm)",'Interpreter','latex')
ylabel('Probability density','Interpreter','latex')
fontsize(gcf,scale=1.8)

%% plot PSF image

cfg = 1;
ch = 4;

[allmufit,~] = get_PoissonRate(p,ThetaStore(:,cfg));

T1 = squeeze(allmu(:,:,1,:,cfg));
T2 = squeeze(allspots(:,:,1,:,cfg));
T3 = squeeze(allmufit(:,:,1,:));
%%%% Sum of all excitation channels
% T1 = squeeze(allmu(:,:,1,1:4,cfg)+allmu(:,:,1,5:8,cfg)+allmu(:,:,1,9:12,cfg));
% T2 = squeeze(allspots(:,:,1,1:4,cfg)+allspots(:,:,1,5:8,cfg)+allspots(:,:,1,9:12,cfg));
% T3 = squeeze(allmufit(:,:,1,1:4)+allmufit(:,:,1,5:8)+allmufit(:,:,1,9:12));
%%%% Sum of 4 channels
% T1 = squeeze(allmu(:,:,1,[1 5 9],cfg)+allmu(:,:,1,[2 6 10],cfg)+allmu(:,:,1,[3 7 11],cfg)+allmu(:,:,1,[4 8 12],cfg));
% T2 = squeeze(allspots(:,:,1,[1 5 9],cfg)+allspots(:,:,1,[2 6 10],cfg)+allspots(:,:,1,[3 7 11],cfg)+allspots(:,:,1,[4 8 12],cfg));
% T3 = squeeze(allmufit(:,:,1,[1 5 9])+allmufit(:,:,1,[2 6 10])+allmufit(:,:,1,[3 7 11])+allmufit(:,:,1,[4 8 12]));

cmax = ceil(max([max(T1(:)),max(T2(:)),max(T3(:))]));
%%
figure
imagesc([-1,1],[-1,1],T1(:,:,ch));
set(gca,'ydir','normal');
colormap("hot");clim([0,cmax])
axis equal
axis off
copygraphics(gca,'ContentType','vector')
%%
figure
imagesc([-1,1],[-1,1],T2(:,:,ch));
set(gca,'ydir','normal');    
colormap("hot");clim([0,cmax])
axis equal
axis off
copygraphics(gca,'ContentType','vector')
%%
figure
imagesc([-1,1],[-1,1],T3(:,:,ch));
set(gca,'ydir','normal');
colormap("hot");clim([0,cmax])
axis equal
axis off
copygraphics(gca,'ContentType','vector')
%% colorbar
figure
colormap("hot");clim([0,cmax])
c = colorbar;
set(gca,"Visible","off")
copygraphics(gca,'ContentType','vector')
c.Location = 'south';
c.FontSize = 16;
c.FontWeight = "bold";
c.AxisLocation = "out";
copygraphics(gca,'ContentType','vector')