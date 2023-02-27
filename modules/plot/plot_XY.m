% plot range-angle heatmap
function [axh] = plot_XY(Xcube,rng_grid,agl_grid,cnt)
Nr = size(Xcube,1);   %%%length of Chirp(num of rangeffts)
Ne = size(Xcube,2);   %%%number of angleffts
Nd = size(Xcube,3);   %%%length of chirp loop

Xpow = abs(Xcube);
Xpow = squeeze(sum(Xpow,3)/size(Xpow,3));

% noisefloor = db2pow(-15);
Xsnr=Xpow;

x_grid=-4:0.03:1;
y_grid=0:0.03:16;

Xsnr_xy=zeros(length(y_grid),length(x_grid));

for i=1:length(x_grid)
    for j=1:length(y_grid)
        r=sqrt(x_grid(i)^2+y_grid(j)^2);
        theta=atan(x_grid(i)/y_grid(j))*180/pi;
        
        if r>rng_grid(end)
            Xsnr_xy(j,i)=0;
        else
            [~,r_idx]=min(abs(rng_grid-r));
            [~,theta_idx]=min(abs(agl_grid-theta));
            Xsnr_xy(j,i)=Xsnr(r_idx,theta_idx);
        end
    end
end

%% plot 2D(range-angle)

% Xsnr = pow2db(Xpow/noisefloor);

% figure('visible','on')
% figure()
% set(gcf,'Position',[10,10,530,420])
[axh] = surf(x_grid,y_grid,Xsnr_xy);
view(0,90)

colorbar
% caxis([0 50])
axis([-4 1 0 16]);
grid off
shading interp

xlabel('X(meters)')
ylabel('Y(meters)')

% caxis([0,5*10^5])
title(strcat('X-Y heatmap, Frame no = ',num2str(cnt)))


end