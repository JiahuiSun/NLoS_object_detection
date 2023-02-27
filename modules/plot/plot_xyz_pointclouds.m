% plot 3D point clouds
function [axh] = plot_xyz_pointclouds(detout,static_idx,dynamic_idx,cnt,Pfa)
% detout format: % [range bin, velocity bin, angle bin, power, range(m), ...
% velocity (m/s), angle(degree)]

% figure('visible','on')
% x-direction: Doppler, y-direction: angle, z-direction: 

% detout_sorted = sortrows(detout,4,'descend');

power = detout(:,4);
x_value=detout(:,8);
y_value=detout(:,9);
z_value=detout(:,10);

% [axh] = scatter3(x_value,y_value,z_value,20,'b', 'filled');
[axh] = scatter3(x_value(static_idx),y_value(static_idx),z_value(static_idx),20,'b', 'filled');
hold on
scatter3(x_value(dynamic_idx),y_value(dynamic_idx),z_value(dynamic_idx),20,'r', 'filled');
hold off

% [axh] = scatter3(x_value(1:ceil(size(detout_sorted,1)/2)),y_value(1:ceil(size(detout_sorted,1)/2)),z_value(1:ceil(size(detout_sorted,1)/2)),20,'r', 'filled');
% hold on
% scatter3(x_value(ceil(size(detout_sorted,1)/2):size(detout_sorted,1)),y_value(ceil(size(detout_sorted,1)/2):size(detout_sorted,1)),z_value(ceil(size(detout_sorted,1)/2):size(detout_sorted,1)),20,'b', 'filled');
% hold off

ax = gca;
% ax.XDir = 'reverse';
view(0,90)
% cb = colorbar;                                     
% cb.Label.String = 'power';


xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis([-5 5 0 10 -5 5]);
title(strcat('xyz point clouds, Frame no = ',num2str(cnt),', Pfa = ',num2str(Pfa)))
grid on

end