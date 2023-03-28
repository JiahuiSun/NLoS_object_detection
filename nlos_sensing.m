% 根据雷达自身位置和转角参数，计算映射前的NLoS区域的理论上界，用该上界过滤点云，并对剩余的点映射
% input: 
%   point_cloud: [x, y, vel] Nx3
%   radar_pos: [x, y]
%   corner_args: struct, parameters for the corner
% output: point_cloud_filter
function point_cloud_filter = nlos_sensing(point_cloud, radar_pos, corner_args)
    point_cloud_ext = [point_cloud(:, 1:2), ones([size(point_cloud, 1), 1])];
    top_wall_y = corner_args.top_wall_y;
    inner_corner = corner_args.inner_corner;
    % 用前墙做NLoS感知
    top_map_corner = [inner_corner(1), 2*top_wall_y-inner_corner(2)];
    top_map_radar = [radar_pos(1), 2*top_wall_y-radar_pos(2)];
    top_border = top_map_corner(2);
    if inner_corner(1) < radar_pos(1)
        left_border = line_by_2p(radar_pos, inner_corner);
        right_border = line_by_2p(top_map_radar, top_map_corner);
    elseif inner_corner(1) > radar_pos(1)
        left_border = line_by_2p(top_map_radar, top_map_corner);
        right_border = line_by_2p(radar_pos, inner_corner);
    end
    flag1 = point_cloud_ext*left_border > 0;
    flag2 = point_cloud_ext*right_border < 0;
    flag3 = point_cloud(:, 2) < top_border;
    flag = bitand(bitand(flag1, flag2), flag3);
    % 映射
    point_cloud(flag, 2) = 2*top_wall_y - point_cloud(flag, 2);
    point_cloud_filter = point_cloud(flag, :);
end
