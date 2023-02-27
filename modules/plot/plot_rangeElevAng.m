function [axh] = plot_rangeElevAng(elevAngdata,rng_grid,vel_grid)

figure('visible','on')
imagesc([rng_grid(1) rng_grid(end)],[vel_grid(1) vel_grid(end)],elevAngdata)
colorbar
caxis([-180 180])
axis([0 5 -1 1]);
set(gca,'YDir','normal')%对Y方向反转
xlabel('Range(meters)')
ylabel('Velocity(m/s)')
title('Elevation Angle')
end