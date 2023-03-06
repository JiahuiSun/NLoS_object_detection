% radar_xy: 雷达坐标系下的点坐标，Nx2
% delta_x, delta_y: 雷达在世界坐标系中的坐标
% yaw: 逆时针从雷达坐标系旋转多少度与世界坐标系重合
function [world_xy] = transform(radar_xy, delta_x, delta_y, yaw)
yaw = yaw * pi / 180;
rotation_matrix = [cos(yaw), sin(yaw); -sin(yaw), cos(yaw)];
translation_vector = [delta_x; delta_y];
% 创建12x2的矩阵并赋值
translation_matrix = repmat(translation_vector, 1, size(radar_xy, 1));
world_xy = (rotation_matrix*(radar_xy') + translation_matrix)';
end