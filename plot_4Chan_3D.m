function plot_4Chan_3D(p,PSF)

PSF = PSF./max(PSF(:));

channel_name = ["x1","y1","x2","y2"];
% channel_name = ["0°","90°","45°","135°"];
% channel_name = ["x1(L)","x2(L)","y1(R)","y2(R)"];

figure
set(gcf,'position',[100 100 1100 700])
ha = tight_subplot(3,5,[0.05 0.04],[0.06 0.1],0.04);
% sgtitle(f_title,'FontName','Arial','FontSize',16,'FontWeight','bold')

Lm = p.xl(end)*1e6;
LmI = 0.5;  % um, size of cross sections

Zm = p.Lz/2*1e6;

for ii=1:4
    h1=ii; h2=ii+5; h3=ii+10;
    axes(ha(h1))
    hold on
    imagesc([-Lm,Lm],[-Lm,Lm],PSF(:,:,ceil(p.Nz/2),ii))
    set(gca,'YDir','normal')
    clim([0,1])
    plot([0 0],[-1 1],'-.r','LineWidth',1.2)
    plot([-1 1],[0 0],'-.r','LineWidth',1.2)
    set(gca,'FontWeight','bold')
    xlabel('y (µm)')
    ylabel('x (µm)')
    % axis equal
    axis([-LmI,LmI,-LmI,LmI])
    ha(h1).XTick = [-LmI,0,LmI]; ha(h1).YTick = [-LmI,0,LmI];

    axes(ha(h2))
    hold on
    imagesc([-Lm,Lm],[-Zm,Zm],squeeze(PSF(:,ceil(p.Nx/2),:,ii))')
    set(gca,'YDir','normal')
    clim([0,1])
    plot([0 0],[-1 1],'-.r','LineWidth',1.2)
    plot([-1 1],[0 0],'-.r','LineWidth',1.2)
    set(gca,'FontWeight','bold')
    xlabel('x (µm)')
    ylabel('z (µm)')
    % axis equal
    axis([-LmI,LmI,-Zm,Zm])
    ha(h2).XTick = [-LmI,0,LmI]; ha(h2).YTick = [-Zm,0,Zm];

    axes(ha(h3))
    hold on
    imagesc([-Lm,Lm],[-Zm,Zm],squeeze(PSF(ceil(p.Nx/2),:,:,ii))')
    set(gca,'YDir','normal')
    clim([0,1])
    plot([0 0],[-1 1],'-.r','LineWidth',1.2)
    plot([-1 1],[0 0],'-.r','LineWidth',1.2)
    set(gca,'FontWeight','bold')
    xlabel('y (µm)')
    ylabel('z (µm)')
    % axis equal
    axis([-LmI,LmI,-Zm,Zm])
    ha(h3).XTick = [-LmI,0,LmI]; ha(h3).YTick = [-Zm,0,Zm];
end

PSF_y = squeeze(PSF(ceil(p.Nx/2),:,ceil(p.Nz/2),:));
PSF_x = squeeze(PSF(:,ceil(p.Nx/2),ceil(p.Nz/2),:));
PSF_z = squeeze(PSF(ceil(p.Nx/2),ceil(p.Nx/2),:,:));
zlm = p.zl*1e6;
xlm = p.xl*1e6;

axes(ha(1*5))
plot(zlm,PSF_z(:,1),'b',zlm,PSF_z(:,2),'r',zlm,PSF_z(:,3),'-.g',zlm,PSF_z(:,4),'--m','LineWidth',1.2)
xlim([zlm(1),zlm(end)])
ylim([0,1])
grid
legend(channel_name(1),channel_name(2),channel_name(3),channel_name(4))
legend('boxoff')
set(gca,'FontWeight','bold')
xlabel('z (µm)')

axes(ha(2*5))
plot(xlm,PSF_x(:,1),'b',xlm,PSF_x(:,2),'r',xlm,PSF_x(:,3),'-.g',xlm,PSF_x(:,4),'--m','LineWidth',1.2)
xlim([-LmI,LmI])
ylim([0,1])
grid
legend(channel_name(1),channel_name(2),channel_name(3),channel_name(4))
legend('boxoff')
set(gca,'FontWeight','bold')
xlabel('x (µm)')

axes(ha(3*5))
plot(xlm,PSF_y(:,1),'b',xlm,PSF_y(:,2),'r',xlm,PSF_y(:,3),'-.g',xlm,PSF_y(:,4),'--m','LineWidth',1.2)
xlim([-LmI,LmI])
ylim([0,1])
grid
legend(channel_name(1),channel_name(2),channel_name(3),channel_name(4))
legend('boxoff')
set(gca,'FontWeight','bold')
xlabel('y (µm)')

for ii=1:4
title(ha(ii),channel_name(ii),'FontName','Arial','FontSize',14,'FontWeight','bold'); 
end
title(ha(5),'Normalized intensity','FontName','Arial','FontSize',14,'FontWeight','bold'); 



