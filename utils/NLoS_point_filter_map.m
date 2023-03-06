% 根据雷达自身位置和转角参数，计算映射前的NLoS区域的理论上界，用该上界过滤点云，并对剩余的点映射
% input: 
%   point_cloud: [x, y, vel] Nx3
%   radar_pos: [x, y]
%   corner_type: 1-12
%   corner_args: struct, parameters for the corner
% output: point_cloud_filter
function point_cloud_filter = NLoS_point_filter_map(point_cloud, radar_pos, corner_type, corner_args)
    point_cloud_ext = [point_cloud(:, 1:2), ones([size(point_cloud, 1), 1])];
    switch corner_type
        case 1
            top_wall_y = corner_args.top_wall_y;
            bottom_wall_y = corner_args.bottom_wall_y;
            left_wall_x = corner_args.left_wall_x;
            right_wall_x = corner_args.right_wall_x;
            
            % 用前墙做NLoS感知
            right_map_radar = [2*right_wall_x-radar_pos(1), radar_pos(2)];
            right_map_left_x = 2*right_wall_x - left_wall_x;
            top_border1 = line_by_2p(right_map_radar, [right_map_left_x, bottom_wall_y]);
            bottom_border1 = line_by_2p(radar_pos, [left_wall_x, bottom_wall_y]);
            right_border1 = right_map_left_x;
            flag11 = point_cloud_ext*bottom_border1 > 0;
            flag12 = point_cloud_ext*top_border1 > 0;
            flag13 = point_cloud_ext < right_border1;
            flag1 = bitand(bitand(flag11, flag12), flag13);
            % 映射
            point_cloud(flag1, 1) = 2*right_wall_x - point_cloud(flag1, 1);

            % 用侧墙做NLoS感知
            top_map_radar = [radar_pos(1), 2*top_wall_y-radar_pos(2)];
            top_map_bottom_y = 2*top_wall_y - bottom_wall_y;
            top_border2 = line_by_2p(radar_pos, [left_wall_x, top_map_bottom_y]);
            bottom_border2 = line_by_2p(top_map_radar, [left_wall_x, top_map_bottom_y]);
            right_border2 = right_wall_x;
            flag21 = point_cloud_ext*top_border2 > 0;
            flag22 = point_cloud_ext*bottom_border2 < 0;
            flag23 = point_cloud_ext(:, 1) < right_border2;
            flag2 = bitand(bitand(flag21, flag22), flag23);
            % 映射
            point_cloud(flag2, 2) = 2*top_wall_y - point_cloud(flag2, 2);

            flag = bitor(flag1, flag2);
        case 2
            top_wall_y = corner_args.top_wall_y;
            bottom_wall_y = corner_args.bottom_wall_y;
            left_wall_x = corner_args.left_wall_x;
            right_wall_x = corner_args.right_wall_x;

            % 用前墙做NLoS感知
            top_map_radar = [radar_pos(1), 2*top_wall_y-radar_pos(2)];
            top_map_bottom_y = 2*top_wall_y - bottom_wall_y;
            left_border1 = line_by_2p(radar_pos, [left_wall_x, bottom_wall_y]);
            right_border1 = line_by_2p(top_map_radar, [left_wall_x, top_map_bottom_y]);
            top_border1 = top_map_bottom_y;
            flag11 = point_cloud_ext*left_border1 > 0;
            flag12 = point_cloud_ext*right_border1 < 0;
            flag13 = point_cloud(:, 2) < top_border1;
            flag1 = bitand(bitand(flag11, flag12), flag13);
            % 映射
            point_cloud(flag1, 2) = 2*top_wall_y - point_cloud(flag1, 2);

            % 用侧墙做NLoS感知
            right_map_left_x = 2*right_wall_x - left_wall_x;
            right_map_radar = [2*right_wall_x - radar_pos(1), radar_pos(2)];
            left_border2 = line_by_2p(right_map_radar, [right_map_left_x, bottom_wall_y]);
            right_border2 = line_by_2p(radar_pos, [right_map_left_x, bottom_wall_y]);
            top_border2 = top_wall_y;
            flag21 = point_cloud_ext*left_border2 > 0;
            flag22 = point_cloud_ext*right_border2 < 0;
            flag23 = point_cloud(:, 2) < top_border2;
            flag2 = bitand(bitand(flag21, flag22), flag23);
            % 映射
            point_cloud(flag2, 1) = 2*right_wall_x - point_cloud(flag2, 1);

            flag = bitor(flag1, flag2);
        case 3
            top_wall_y = corner_args.top_wall_y;
            bottom_wall_y = corner_args.bottom_wall_y;
            left_wall_x = corner_args.left_wall_x;

            % 用侧墙做NLoS感知
            top_map_bottom_y = 2*top_wall_y - bottom_wall_y;
            top_map_radar = [radar_pos(1), 2*top_wall_y-radar_pos(2)];
            left_border = line_by_2p(radar_pos, [left_wall_x, top_map_bottom_y]);
            right_border = line_by_2p(top_map_radar, [left_wall_x, top_map_bottom_y]);
            flag1 = point_cloud_ext*left_border > 0;
            flag2 = point_cloud_ext*right_border < 0;
            flag = bitand(flag1, flag2);
            % 映射
            point_cloud(flag, 2) = 2*top_wall_y - point_cloud(flag, 2);
        case 4
            top_wall_y = corner_args.top_wall_y;
            bottom_wall_y = corner_args.bottom_wall_y;
            left_wall_x = corner_args.left_wall_x;

            % 用前墙做NLoS感知
            top_map_bottom_y = 2*top_wall_y - bottom_wall_y;
            top_map_radar = [radar_pos(1), 2*top_wall_y-radar_pos(2)];
            left_border = line_by_2p(radar_pos, [left_wall_x, bottom_wall_y]);
            right_border = line_by_2p(top_map_radar, [left_wall_x, top_map_bottom_y]);
            top_border = top_map_bottom_y;
            flag1 = point_cloud_ext*left_border > 0;
            flag2 = point_cloud_ext*right_border < 0;
            flag3 = point_cloud(:, 2) < top_border;
            flag = bitand(bitand(flag1, flag2), flag3);
            % 映射
            point_cloud(flag, 2) = 2*top_wall_y - point_cloud(flag, 2);
        case 5
            top_wall_y = corner_args.top_wall_y;
            bottom_wall_y = corner_args.bottom_wall_y;
            right_wall_x = corner_args.left_wall_x;

            % 用侧墙做NLoS感知
            top_map_radar = [radar_pos(1), 2*top_wall_y-radar_pos(2)];
            top_map_bottom_y = 2*top_wall_y - bottom_wall_y;
            left_border = line_by_2p(top_map_radar, [right_wall_x, top_map_bottom_y]);
            right_border = line_by_2p(radar_pos, [right_wall_x, top_map_bottom_y]);
            flag1 = point_cloud_ext*left_border > 0;
            flag2 = point_cloud_ext*right_border < 0;
            flag = bitand(flag1, flag2);
            % 映射
            point_cloud(flag, 2) = 2*top_wall_y - point_cloud(flag, 2);
        case 6
            top_wall_y = corner_args.top_wall_y;
            bottom_wall_y = corner_args.bottom_wall_y;
            right_wall_x = corner_args.left_wall_x;

            % 用前墙做NLoS感知
            top_map_radar = [radar_pos(1), 2*top_wall_y-radar_pos(2)];
            top_map_bottom_y = 2*top_wall_y - bottom_wall_y;
            left_border = line_by_2p(top_map_radar, [right_wall_x, top_map_bottom_y]);
            right_border = line_by_2p(radar_pos, [right_wall_x, bottom_wall_y]);
            top_border = top_map_bottom_y;
            flag1 = point_cloud_ext*left_border > 0;
            flag2 = point_cloud_ext*right_border < 0;
            flag3 = point_cloud_ext(:, 2) < top_border;
            flag = bitand(bitand(flag1, flag2), flag3);
            % 映射
            point_cloud(flag, 2) = 2*top_wall_y - point_cloud(flag, 2);
        case 7
            % NOTE: 和1一模一样
            top_wall_y = corner_args.top_wall_y;
            bottom_wall_y = corner_args.bottom_wall_y;
            left_wall_x = corner_args.left_wall_x;
            right_wall_x = corner_args.right_wall_x;

            % 用前墙
            right_map_radar = [2*right_wall_x-radar_pos(1), radar_pos(2)];
            right_map_left_x = 2*right_wall_x - left_wall_x;
            top_border1 = line_by_2p(right_map_radar, [right_map_left_x, bottom_wall_y]);
            bottom_border1 = line_by_2p(radar_pos, [left_wall_x, bottom_wall_y]);
            right_border1 = right_map_left_x;
            flag11 = point_cloud_ext*bottom_border1 > 0;
            flag12 = point_cloud_ext*top_border1 > 0;
            flag13 = point_cloud_ext < right_border1;
            flag1 = bitand(bitand(flag11, flag12), flag13);
            % 映射
            point_cloud(flag1, 1) = 2*right_wall_x - point_cloud(flag1, 1);

            % 用侧墙
            top_map_radar = [radar_pos(1), 2*top_wall_y-radar_pos(2)];
            top_map_bottom_y = 2*top_wall_y - bottom_wall_y;
            top_border2 = line_by_2p(radar_pos, [left_wall_x, top_map_bottom_y]);
            bottom_border2 = line_by_2p(top_map_radar, [left_wall_x, top_map_bottom_y]);
            right_border2 = right_wall_x;
            flag21 = point_cloud_ext*top_border2 > 0;
            flag22 = point_cloud_ext*bottom_border2 < 0;
            flag23 = point_cloud_ext(:, 1) < right_border2;
            flag2 = bitand(bitand(flag21, flag22), flag23);
            % 映射
            point_cloud(flag2, 2) = 2*top_wall_y - point_cloud(flag2, 2);

            flag = bitor(flag1, flag2);
        case 8
            % NOTE: 和1完全对称
            top_wall_y = corner_args.top_wall_y;
            bottom_wall_y = corner_args.bottom_wall_y;
            left_wall_x = corner_args.left_wall_x;
            right_wall_x = corner_args.right_wall_x;

            % 用前墙
            left_map_radar = [2*left_wall_x-radar_pos(1), radar_pos(2)];
            left_map_right_x = 2*left_wall_x - right_wall_x;
            top_border1 = line_by_2p(left_map_radar, [left_map_right_x, bottom_wall_y]);
            bottom_border1 = line_by_2p(radar_pos, [right_wall_x, bottom_wall_y]);
            left_border1 = left_map_right_x;
            flag11 = point_cloud_ext*bottom_border1 < 0;
            flag12 = point_cloud_ext*top_border1 < 0;
            flag13 = point_cloud_ext > left_border1;
            flag1 = bitand(bitand(flag11, flag12), flag13);
            % 映射
            point_cloud(flag1, 1) = 2*left_wall_x - point_cloud(flag1, 1);

            % 用侧墙
            top_map_radar = [radar_pos(1), 2*top_wall_y-radar_pos(2)];
            top_map_bottom_y = 2*top_wall_y - bottom_wall_y;
            top_border2 = line_by_2p(radar_pos, [right_wall_x, top_map_bottom_y]);
            bottom_border2 = line_by_2p(top_map_radar, [right_wall_x, top_map_bottom_y]);
            left_border2 = left_wall_x;
            flag21 = point_cloud_ext*top_border2 < 0;
            flag22 = point_cloud_ext*bottom_border2 > 0;
            flag23 = point_cloud_ext(:, 1) > left_border2;
            flag2 = bitand(bitand(flag21, flag22), flag23);
            % 映射
            point_cloud(flag2, 2) = 2*top_wall_y - point_cloud(flag2, 2);

            flag = bitor(flag1, flag2);
        case 9
            top_wall_y = corner_args.top_wall_y;
            bottom_wall_y = corner_args.bottom_wall_y;
            left_wall_x = corner_args.left_wall_x;
            right_wall_x = corner_args.right_wall_x;

            % NOTE: 和2一模一样
            % 用前墙看左边
            top_map_radar = [radar_pos(1), 2*top_wall_y-radar_pos(2)];
            top_map_bottom_y = 2*top_wall_y - bottom_wall_y;
            left_border1 = line_by_2p(radar_pos, [left_wall_x, bottom_wall_y]);
            right_border1 = line_by_2p(top_map_radar, [left_wall_x, top_map_bottom_y]);
            top_border1 = top_map_bottom_y;
            flag11 = point_cloud_ext*left_border1 > 0;
            flag12 = point_cloud_ext*right_border1 < 0;
            flag13 = point_cloud(:, 2) < top_border1;
            flag1 = bitand(bitand(flag11, flag12), flag13);
            % 映射
            point_cloud(flag1, 2) = 2*top_wall_y - point_cloud(flag1, 2);

            % 用右墙看左边
            right_map_left_x = 2*right_wall_x - left_wall_x;
            right_map_radar = [2*right_wall_x - radar_pos(1), radar_pos(2)];
            left_border2 = line_by_2p(right_map_radar, [right_map_left_x, bottom_wall_y]);
            right_border2 = line_by_2p(radar_pos, [right_map_left_x, bottom_wall_y]);
            top_border2 = top_wall_y;
            flag21 = point_cloud_ext*left_border2 > 0;
            flag22 = point_cloud_ext*right_border2 < 0;
            flag23 = point_cloud(:, 2) < top_border2;
            flag2 = bitand(bitand(flag21, flag22), flag23);
            % 映射
            point_cloud(flag2, 1) = 2*right_wall_x - point_cloud(flag2, 1);

            % NOTE: 和2对称，和6一模一样
            % 用前墙看右边
            top_map_radar = [radar_pos(1), 2*top_wall_y-radar_pos(2)];
            top_map_bottom_y = 2*top_wall_y - bottom_wall_y;
            left_border = line_by_2p(top_map_radar, [right_wall_x, top_map_bottom_y]);
            right_border = line_by_2p(radar_pos, [right_wall_x, bottom_wall_y]);
            top_border = top_map_bottom_y;
            flag31 = point_cloud_ext*left_border > 0;
            flag32 = point_cloud_ext*right_border < 0;
            flag33 = point_cloud_ext(:, 2) < top_border;
            flag3 = bitand(bitand(flag31, flag32), flag33);
            % 映射
            point_cloud(flag3, 2) = 2*top_wall_y - point_cloud(flag3, 2);

            % 用左墙看右边
            left_map_right_x = 2*left_wall_x - right_wall_x;
            left_map_radar = [2*left_wall_x - radar_pos(1), radar_pos(2)];
            left_border2 = line_by_2p(radar_pos, [left_map_right_x, bottom_wall_y]);
            right_border2 = line_by_2p(left_map_radar, [left_map_right_x, bottom_wall_y]);
            top_border2 = top_wall_y;
            flag41 = point_cloud_ext*left_border2 > 0;
            flag42 = point_cloud_ext*right_border2 < 0;
            flag43 = point_cloud(:, 2) < top_border2;
            flag4 = bitand(bitand(flag41, flag42), flag43);
            % 映射
            point_cloud(flag4, 1) = 2*left_wall_x - point_cloud(flag4, 1);

            flag = bitor(bitor(bitor(flag1, flag2), flag3), flag4);
        case 10
            % NOTE: 和7一模一样
            bottom_wall_y = corner_args.bottom_wall_y;
            left_wall_x = corner_args.left_wall_x;
            right_wall_x = corner_args.right_wall_x;

            % 用前墙
            right_map_radar = [2*right_wall_x-radar_pos(1), radar_pos(2)];
            right_map_left_x = 2*right_wall_x - left_wall_x;
            top_border1 = line_by_2p(right_map_radar, [right_map_left_x, bottom_wall_y]);
            bottom_border1 = line_by_2p(radar_pos, [left_wall_x, bottom_wall_y]);
            right_border1 = right_map_left_x;
            flag1 = point_cloud_ext*bottom_border1 > 0;
            flag2 = point_cloud_ext*top_border1 > 0;
            flag3 = point_cloud_ext < right_border1;
            flag = bitand(bitand(flag1, flag2), flag3);
            % 映射
            point_cloud(flag, 1) = 2*right_wall_x - point_cloud(flag, 1);
        case 11
            % NOTE: 和8一模一样
            bottom_wall_y = corner_args.bottom_wall_y;
            left_wall_x = corner_args.left_wall_x;
            right_wall_x = corner_args.right_wall_x;

            % 用前墙
            left_map_radar = [2*left_wall_x-radar_pos(1), radar_pos(2)];
            left_map_right_x = 2*left_wall_x - right_wall_x;
            top_border1 = line_by_2p(left_map_radar, [left_map_right_x, bottom_wall_y]);
            bottom_border1 = line_by_2p(radar_pos, [right_wall_x, bottom_wall_y]);
            left_border1 = left_map_right_x;
            flag1 = point_cloud_ext*bottom_border1 < 0;
            flag2 = point_cloud_ext*top_border1 < 0;
            flag3 = point_cloud_ext > left_border1;
            flag = bitand(bitand(flag1, flag2), flag3);
            % 映射
            point_cloud(flag, 1) = 2*left_wall_x - point_cloud(flag, 1);
        case 12
            % NOTE: 和3和5类似
            bottom_wall_y = corner_args.bottom_wall_y;
            left_wall_x = corner_args.left_wall_x;
            right_wall_x = corner_args.right_wall_x;

            % 用左墙看右边
            left_map_radar = [2*left_wall_x-radar_pos(1), radar_pos(2)];
            left_map_right_x = 2*left_wall_x - right_wall_x;
            left_border1 = line_by_2p(radar_pos, [left_map_right_x, bottom_wall_y]);
            right_border1 = line_by_2p(left_map_radar, [left_map_right_x, bottom_wall_y]);
            flag11 = point_cloud_ext*left_border1 > 0;
            flag12 = point_cloud_ext*right_border1 < 0;
            flag1 = bitand(flag11, flag12);
            % 映射
            point_cloud(flag1, 1) = 2*left_wall_x - point_cloud(flag1, 1);

            % 用右墙看左边
            right_map_left_x = 2*right_wall_x - left_wall_x;
            right_map_radar = [2*right_wall_x-radar_pos(1), radar_pos(2)];
            left_border2 = line_by_2p(right_map_radar, [right_map_left_x, bottom_wall_y]);
            right_border2 = line_by_2p(radar_pos, [right_map_left_x, bottom_wall_y]);
            flag21 = point_cloud_ext*left_border2 > 0;
            flag22 = point_cloud_ext*right_border2 < 0;
            flag2 = bitand(flag21, flag22);
            % 映射
            point_cloud(flag2, 1) = 2*right_wall_x - point_cloud(flag2, 1);

            flag = bitor(flag1, flag2);
        otherwise
            error("Such a corner type is not considered.");
    end
    point_cloud_filter = point_cloud(flag, :);
end
