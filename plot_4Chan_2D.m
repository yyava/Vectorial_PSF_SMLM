function plot_4Chan_2D(p,images)
% Intput: 
%   images(Nx,Nx,(Nz),Nc)

images = squeeze(images);
channel_name = ["x1","y1","x2","y2"];
% channel_name = ["0°","90°","45°","135°"];
% channel_name = ["x1 (L)","y1 (R)","x2 (L)","y2 (R)"];

figure

set(gcf,'position',[100 100 500 500])
ha = tight_subplot(2,2,[0.06 0.02],[0.06 0.12],0.05);

% s_theta = ['\theta = ',num2str(p.pola/pi*180,'%d'),'° '];
% s_phi = ['\phi = ',num2str(p.azim/pi*180,'%d'),'° '];
% orient_angle = [s_theta,s_phi];
% f_title =strcat("Half wave plate z=80 nm, ",orient_angle);
% sgtitle(f_title,'FontName','Arial','FontSize',16,'FontWeight','bold')

Lm = p.xl(end)*1e9;
muMax = max(images(:));
muMin = min(images(:));

for ii=1:4
    axes(ha(ii))
    imagesc([-Lm,Lm],[-Lm,Lm],images(:,:,ii))
    set(gca,'YDir','normal')
    title(channel_name(ii))
    axis('square')
    axis off
    fontsize(gca,scale=1.5)
    % colorbar
    clim([muMin,muMax])

end

end