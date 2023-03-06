clc
clearvars
close all

addpath(genpath('config'));
addpath(genpath('utils'));
addpath(genpath('modules/plot'));
addpath(genpath('modules/fft'));
addpath(genpath('modules/detection'));
addpath(genpath('modules/cluster'));

%% 设置参数
params = get_params_value_2();

% constant parameters
c = params.c; % Speed of light in air (m/s)
fc = params.fc; % Center frequency (Hz)
lambda = params.lambda;
Rx = params.Rx;
Tx = params.Tx;

% configuration parameters
Fs = params.Fs;
sweepSlope = params.sweepSlope;
samples = params.samples;
loop = params.loop;

Tc = params.Tc; % us 
fft_Rang = params.fft_Rang;
fft_Vel = params.fft_Vel;
fft_Ang = params.fft_Ang;
num_crop = params.num_crop;
max_value = params.max_value; % normalization the maximum of data WITH 1843

% Creat grid table
rng_grid = params.rng_grid;
agl_grid = params.agl_grid;
vel_grid = params.vel_grid;

% Algorithm parameters
data_each_frame = samples*loop*Tx;
Is_Windowed = 1;% 1==> Windowing before doing range and angle fft

data_path = "/home/agent/下载/";
radar_pos_key = "1-0";
person_pos_key = "2_0";
radar_angle = 45;
corner_type = 4;
n_frame = 30;

switch person_pos_key
    case "1_0"
        ground_truth = [-3.6, 4.8];
    case "1_1"
        ground_truth = [-4.8, 4.8];
    case "1_2"
        ground_truth = [-6.0, 4.8];
    case "2_0"
        ground_truth = [-3.6, 3.6];
    case "2_1"
        ground_truth = [-4.8, 3.6];
    case "2_2"
        ground_truth = [-6.0, 3.6];
    otherwise
        error("wrong person position")
end

switch radar_pos_key
    case "1-0"
        radar_pos = [-0.6, 2.4];
    case "1-1"
        radar_pos = [-0.6, 1.2];
    case "1-2"
        radar_pos = [-0.6, 0.0];
    case "2-0"
        radar_pos = [0.0, 2.4];
    case "2-1"
        radar_pos = [0.0, 1.2];
    case "2-2"
        radar_pos = [0.0, 0.0];
    otherwise
        error("wrong radar position")
end

top_wall_y = 5.32;
bottom_wall_y = 3.18;
left_wall_x = -1.8;
fov_line_k = (bottom_wall_y - radar_pos(2)) / (left_wall_x - radar_pos(1));
fov_line_z = radar_pos(2) - fov_line_k * radar_pos(1);
switch corner_type
    case 1
    case 2
    case 3
        corner_args = struct('top_wall_y', 5.32, 'bottom_wall_y', 3.18, 'left_wall_x', -1.8);
    case 4
        corner_args = struct('top_wall_y', 5.32, 'bottom_wall_y', 3.18, 'left_wall_x', -1.8);
    case 5
    case 6
    case 7 
    case 8
    case 9
    case 10
    case 11
    case 12
    otherwise
        error("Such a corner type is not considered.");
end

%% 读数据，处理每一帧
tic;
filename = data_path + "radar" + radar_pos_key + "_person" + person_pos_key + "_" + radar_angle;
fid = fopen(filename+"_Raw_0.bin", 'r');
all_frame_final_data = [];
all_frame_original_data = [];
final_gif = filename + '.gif';
final_png = filename + ".png";
bar = waitbar(0, "please wait");
for cnt = 1:n_frame
    waitbar(cnt/n_frame, bar, "Running  "+num2str(cnt/n_frame*100, "%.1f")+"%");
    data_frame = get_next_frame(fid, params);
    data_chirp = [];
    for cj = 1:Tx*loop
        temp_data = data_frame(:, (cj-1)*samples+1:cj*samples);
        data_chirp(:, :, cj) = temp_data;
    end
    
    % separate the chirps for TDM-MIMO with 3 TXs
    chirp_tx1 = data_chirp(:, :, 1:3:end);
    chirp_tx3 = data_chirp(:, :, 2:3:end);
    chirp_tx2 = data_chirp(:, :, 3:3:end);

    % permutation with the format [samples, Rx, chirp]
    chirp_tx1 = permute(chirp_tx1, [2,1,3]);
    chirp_tx2 = permute(chirp_tx2, [2,1,3]);
    chirp_tx3 = permute(chirp_tx3, [2,1,3]);

    %% Range FFT for chirps from 3txs
    [Rangedata_tx1] = fft_range(chirp_tx1, fft_Rang, Is_Windowed);
    [Rangedata_tx2] = fft_range(chirp_tx2, fft_Rang, Is_Windowed);
    [Rangedata_tx3] = fft_range(chirp_tx3, fft_Rang, Is_Windowed);

    %% Doppler FFT
    Dopplerdata_tx1 = fft_doppler(Rangedata_tx1, fft_Vel, 0);
    Dopplerdata_tx2 = fft_doppler(Rangedata_tx2, fft_Vel, 0);
    Dopplerdata_tx3 = fft_doppler(Rangedata_tx3, fft_Vel, 0);
    Dopdata_sum = squeeze(mean(abs(Dopplerdata_tx1), 2));  % 对不同接收天线的多普勒结果取平均
    
    %% CFAR detector on Range-Velocity to detect targets 
    % Output format: [doppler index; range index(start from index 1); cell power]
    Pfa = 1e-7; % probability of false alarm
    [Resl_indx] = cfar_RV(Dopdata_sum, fft_Rang, num_crop, Pfa);  % 这里就是找bin里面最高的峰
    detout = peakGrouping(Resl_indx);
    
    % 画图
    t = tiledlayout('flow');
    % Plot range-Doppler image
    nexttile;
    hold on;
    [axh] = surf(vel_grid, rng_grid, Dopdata_sum);
    view(0, 90);
    axis([-2 2 0 15]);
    grid off;
    shading interp;
    xlabel('Doppler Velocity (m/s)');
    ylabel('Range(meters)');
    colorbar;
    title('Range-Doppler heatmap');
    scatter(vel_grid(detout(1, :)), rng_grid(detout(2, :)), 10, 'r', 'filled');
    hold off;

    %% doppler compensation on Rangedata_even using the max-intensity peak on each range bin
    for ri = num_crop+1:fft_Rang-num_crop
        find_idx = find(detout(2, :) == ri);
        if isempty(find_idx)
            continue
        else
            % pick the first larger velocity
            pick_idx = find_idx(1);
            % phase compensation for virtual elements
            pha_comp_term_tx3 = exp(-1i * 2*pi * (detout(1,pick_idx)-fft_Vel/2-1) / (fft_Vel*3));
            Dopplerdata_tx3(ri, :, :) = Dopplerdata_tx3(ri, :, :) * pha_comp_term_tx3;
            pha_comp_term_tx2 = exp(-1i * 4*pi * (detout(1,pick_idx)-fft_Vel/2-1) / (fft_Vel*3));
            Dopplerdata_tx2(ri, :, :) = Dopplerdata_tx2(ri, :, :) * pha_comp_term_tx2;
        end
    end
    Dopdata_merge = [Dopplerdata_tx1, Dopplerdata_tx3, Dopplerdata_tx2];
    
    %% 估计角度，生成点云图
    [x_vector, y_vector, z_vector] = music_xyz(Dopdata_merge, params, detout, Tx, Rx, -60, 60, 1, -15, 15, 1, 1);
    
    % Transform bin index to range/velocity/angle
    Resel_vel = vel_grid(detout(1, :), 1);
    Resel_rng = rng_grid(detout(2, :), 1);
    % 因为距离要向xyz轴投影得到坐标，music输出是投影的系数
    x_vec = -x_vector'.*Resel_rng;
    y_vec = y_vector'.*Resel_rng;
    z_vec = -z_vector'.*Resel_rng;

    %% 将雷达坐标系下的点云转换到世界坐标系下
    world_xy = transform([x_vec, y_vec], radar_pos(1), radar_pos(2), 360-radar_angle);
    world_xyv = [world_xy, Resel_vel];
    all_frame_original_data = [all_frame_original_data; world_xyv];

    % 画原始点云图
    nexttile;
    hold on;
    scatter(world_xy(Resel_vel==0, 1), world_xy(Resel_vel==0, 2), 20, 'b', 'filled');
    scatter(world_xy(Resel_vel~=0, 1), world_xy(Resel_vel~=0, 2), 20, 'r', 'filled');
    xlabel('x (m)');
    ylabel('y (m)');
    title(strcat('Frame no = ', num2str(cnt), ', Pfa = ', num2str(Pfa)));
    line([-20, 5], [top_wall_y, top_wall_y]);
    line([-20, left_wall_x], [bottom_wall_y, bottom_wall_y]);
    line([left_wall_x, left_wall_x], [bottom_wall_y, 0]);
    line([radar_pos(1), (top_wall_y - fov_line_z) / fov_line_k], [radar_pos(2), top_wall_y]);
    plot(ground_truth(1), ground_truth(2), '*g', "MarkerSize", 10);
    plot(radar_pos(1), radar_pos(2), 'dc', "MarkerSize", 10);
    axis([-20 5 0 20]);
    hold off;

    %% 过滤并映射
    world_xyv_nlos = NLoS_point_filter_map(world_xyv, radar_pos, corner_type, corner_args);

    % 画过滤后的点云图
    nexttile;
    hold on;
    if size(world_xyv_nlos, 1) < 1
        fprintf("no point left.\n");
    else
        all_frame_final_data = [all_frame_final_data; world_xyv_nlos];
        scatter(world_xyv_nlos(world_xyv_nlos(:, 3)==0, 1), world_xyv_nlos(world_xyv_nlos(:, 3)==0, 2), 20, 'b', 'filled');
        scatter(world_xyv_nlos(world_xyv_nlos(:, 3)~=0, 1), world_xyv_nlos(world_xyv_nlos(:, 3)~=0, 2), 20, 'r', 'filled');
    end
    xlabel('x (m)');
    ylabel('y (m)');
    title(strcat('Frame no = ', num2str(cnt), ', Pfa = ', num2str(Pfa)));
    line([-20, 5], [top_wall_y, top_wall_y]);
    line([-20, left_wall_x], [bottom_wall_y, bottom_wall_y]);
    line([left_wall_x, left_wall_x], [bottom_wall_y, 0]);
    line([radar_pos(1), (top_wall_y - fov_line_z) / fov_line_k], [radar_pos(2), top_wall_y]);
    plot(ground_truth(1), ground_truth(2), '*g', "MarkerSize", 10);
    plot(radar_pos(1), radar_pos(2), 'dc', "MarkerSize", 10);
    axis([-20 5 0 20]);
    hold off;

    % 叠加多帧点云图
    nexttile;
    hold on;
    if size(all_frame_final_data, 1) > 0
        static_idx = all_frame_final_data(:, 3)==0;
        dynamic_idx = all_frame_final_data(:, 3)~=0;
        scatter(all_frame_final_data(static_idx, 1), all_frame_final_data(static_idx, 2), 20, 'b', 'filled');
        scatter(all_frame_final_data(dynamic_idx, 1), all_frame_final_data(dynamic_idx, 2), 20, 'r', 'filled');
    end
    xlabel('x (m)');
    ylabel('y (m)');
    title("Stacked point cloud");
    line([-20, 5], [top_wall_y, top_wall_y]);
    line([-20, left_wall_x], [bottom_wall_y, bottom_wall_y]);
    line([left_wall_x, left_wall_x], [bottom_wall_y, 0]);
    line([radar_pos(1), (top_wall_y - fov_line_z) / fov_line_k], [radar_pos(2), top_wall_y]);
    plot(ground_truth(1), ground_truth(2), '*g', "MarkerSize", 10);
    plot(radar_pos(1), radar_pos(2), 'dc', "MarkerSize", 10);
    axis([-20 5 0 20]);
    hold off;

    % 保存动图
    [A, map] = rgb2ind(frame2im(getframe(gcf)), 256);
    if cnt == 1
        imwrite(A, map, final_gif, 'gif', 'Loopcount', inf, 'DelayTime', 0.3);
    else
        imwrite(A, map, final_gif, 'gif', 'WriteMode', 'append', 'DelayTime', 0.3);
    end
end
% 保存最后一帧图片
imwrite(frame2im(getframe(gcf)), final_png);
fprintf(' \nIt took %6.2f s. \n', toc);
fclose(fid);
