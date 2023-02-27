% plot 2D point clouds
function [axh] = plot_2d_pointcloud(x_value, y_value, static_idx, dynamic_idx, cnt, Pfa)
[axh] = scatter(x_value(static_idx), y_value(static_idx), 20, 'b', 'filled');
hold on
scatter(x_value(dynamic_idx), y_value(dynamic_idx), 20, 'r', 'filled');
hold off
ax = gca;
xlabel('x (m)')
ylabel('y (m)')
% axis([-5 5 0 10]);
title(strcat('xy point clouds, Frame no = ', num2str(cnt), ', Pfa = ', num2str(Pfa)))
grid on
end