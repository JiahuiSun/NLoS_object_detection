clc
clearvars
close all

addpath(genpath('.\config'));
addpath(genpath('.\utils'));
addpath(genpath('.\modules\plot'));
addpath(genpath('.\modules\fft'));
addpath(genpath('.\modules\detection'));
addpath(genpath('.\modules\cluster'));

%% parameter setting
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
Is_plot_rangeDop = 1;
Is_plot_rangeAng = 0;
radar_stationary=1;

data_path="/home/agent/下载/";
filename=data_path+"radar1_0_person1_0_45_Raw_0.bin";

fid = fopen(filename,'r');

save_det_data_all=[];

tic;
% data_frames=readDCA1000(filename);
% M1=moviein(20);
% filename1='E:\Project\mmWave\Workspace\gifs\cornerRadar\pcmusic_CR4_1p_dynamic_radarmove_45_1_1e-7.gif';
%% 开始处理数据
for cnt=1:100
    disp(cnt)
%     data_frame = data_frames(:, (cnt-1)*data_each_frame+1:cnt*data_each_frame);
    data_frame=get_next_frame(fid,params);
    
%     if cnt<=30
%         continue
%     end
    data_chirp = [];

    for cj = 1:Tx*loop
        temp_data = data_frame(:, (cj-1)*samples+1:cj*samples);
        data_chirp(:,:,cj) = temp_data;
    end
    
    % separate the chirps for TDM-MIMO with 3 TXs
    chirp_tx1 = data_chirp(:,:,1:3:end);%4*1024*128
    chirp_tx3 = data_chirp(:,:,2:3:end);
    chirp_tx2 = data_chirp(:,:,3:3:end);

    % permutation with the format [samples, Rx, chirp]
    chirp_tx1 = permute(chirp_tx1, [2,1,3]);
    chirp_tx2 = permute(chirp_tx2, [2,1,3]);
    chirp_tx3 = permute(chirp_tx3, [2,1,3]);

    % Range FFT for chirps from 3txs
    [Rangedata_tx1] = fft_range(chirp_tx1,fft_Rang,Is_Windowed);
    [Rangedata_tx2] = fft_range(chirp_tx2,fft_Rang,Is_Windowed);
    [Rangedata_tx3] = fft_range(chirp_tx3,fft_Rang,Is_Windowed);

    % Doppler FFT
    Dopplerdata_tx1 = fft_doppler(Rangedata_tx1, fft_Vel, 0);
    Dopplerdata_tx2 = fft_doppler(Rangedata_tx2, fft_Vel, 0);
    Dopplerdata_tx3 = fft_doppler(Rangedata_tx3, fft_Vel, 0);
    Dopdata_sum = squeeze(mean(abs(Dopplerdata_tx1), 2));  % 这里为什么要对不同天线的多普勒结果取平均
    
    % Plot range-Doppler image
    
    if Is_plot_rangeDop
        plot_rangeDop(Dopdata_sum,rng_grid,vel_grid);
    end
    
%     figure()
%     imagesc([vel_grid(1) vel_grid(end)],[rng_grid(1) rng_grid(end)],Dopdata_sum)
%     colorbar
%     axis([-2 2 0 13]);
%     set(gca,'YDir','normal')%对Y方向反转
%     xlabel('Doppler Velocity (m/s)')
%     ylabel('Range(meters)')
    
    
    %% CFAR detector on Range-Velocity to detect targets 
    % Output format: [doppler index, range index(start from index 1), ...
    % cell power]
    Pfa = 1e-7; % probability of false alarm
    [Resl_indx] = cfar_RV(Dopdata_sum, fft_Rang, num_crop, Pfa);  % 这里就是找bin里面最高的峰，也就意味着这里有物体，但可能有多个，需要通过角度fft
    detout = peakGrouping(Resl_indx);  % 这里的原理是什么
    
%     hold on
%     scatter(vel_grid(detout(1,:)),rng_grid(detout(2,:)),10,'r','filled')
    
% 这里是什么操作？为什么要补偿？
    % doppler compensation on Rangedata_even using the max-intensity peak
    % on each range bin
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
    
    Rangedata_merge = [Rangedata_tx1, Rangedata_tx3];  % 为啥只要1、3
    Dopdata_merge=[Dopplerdata_tx1, Dopplerdata_tx3,Dopplerdata_tx2];
    
%     % Angle FFT (Azimuth angle estimation)
%     Angdata = fft_angle(Dopdata_merge(:,1:2 * Rx, :),fft_Ang,Is_Windowed);
%     Angdata_crop = Angdata(num_crop + 1:fft_Rang - num_crop, :, :);
%     [Angdata] = Normalize(Angdata, max_value);
%     [Angdata_crop] = Normalize(Angdata_crop, max_value);
    
    %% point cloud generation，上面做完了range和velocity fft，下面是处理角度
%     [x_vector, y_vector, z_vector]=naive_xyz(Dopdata_merge,detout,Tx,Rx,fft_Ang,Is_Windowed); 
    [x_vector, y_vector, z_vector]=music_xyz(Dopdata_merge,params,detout,Tx,Rx,-60,60,1,-15,15,1,1); % 这里的原理，我也得学明白
    
    % Angle estimation for detected point clouds
    Dopplerdata_merge = permute([Dopplerdata_tx1, Dopplerdata_tx3], [2, 1, 3]);  % 这里才估计角度吗？那上面music_xyz怎么得到xyz的？
    [Resel_agl, ~, rng_excd_list] = angle_estim_dets(detout, Dopplerdata_merge, fft_Vel, ...
        fft_Ang, Rx, 2, num_crop);
    
    % Transform bin index to range/velocity/angle
    Resel_agl_deg = agl_grid(1, Resel_agl)';
    Resel_vel = vel_grid(detout(1,:), 1);
    Resel_rng = rng_grid(detout(2,:), 1);
    
    x_vec=-x_vector'.*Resel_rng;
    y_vec=y_vector'.*Resel_rng;
    z_vec=-z_vector'.*Resel_rng;
    
        
    % save_det data format below
    % [1 range bin, 2 velocity bin, 3 angle bin, 4 power, 5 range(m), 6 velocity (m/s), 7 angle(degree),8 x ,9 y ,10 z]
    save_det_data = [detout(2,:)', detout(1,:)', Resel_agl', detout(3,:)', ...
        Resel_rng, Resel_vel, Resel_agl_deg,x_vec,y_vec,z_vec];
%     disp(save_det_data)

    % delete complex y value
    valid_idx=find(imag(y_vec)==0);
%     save_det_data=save_det_data(valid_idx,:);  人就是静止的，不能删
    fprintf('save_det_data %d\n',size(save_det_data,1));
    
    vel=save_det_data(:,6);
    theta=atan(save_det_data(:,8)./save_det_data(:,9));
%     这是在干嘛？
    if radar_stationary
        static_idx=find(vel==0);
        dynamic_idx=find(vel~=0);
    end
    
%     plot_xyz_pointclouds(save_det_data,static_idx,dynamic_idx,cnt,Pfa);
%     M1(:,end+1)=getframe(gcf);
%     [A,map]=rgb2ind(frame2im(getframe(gcf)),256);
%     if cnt==1
%         imwrite(A,map,filename1,'gif','Loopcount',inf,'DelayTime',0.3);
%     else
%         imwrite(A,map,filename1,'gif','WriteMode','append','DelayTime',0.3);
%     end
    
    
%     idx = dbscan(save_det_data(:,8:10),1,3);
%     idx_set=unique(idx);
%     save_det_data_dbscan=save_det_data(find(idx~=-1),:);
% 
%     save_det_data_all=[save_det_data_all;save_det_data_dbscan];
%     fprintf('save_det_data_dbscan %d\n',size(save_det_data_dbscan,1));
    
    save_det_data_all=[save_det_data_all;save_det_data];
    break
end
fprintf(' \nIt took %6.2f s. \n',toc);

% idx = dbscan(save_det_data_all(:,8:10),1,15);
% idx_set=unique(idx);
% save_det_data_dbscan=save_det_data_all(find(idx~=-1),:);
% 
% % save_det_data_all=[save_det_data_all;save_det_data_dbscan];
% fprintf('save_det_data_dbscan %d\n',size(save_det_data_dbscan,1));

plot_xyz_pointclouds(save_det_data_all,static_idx,dynamic_idx,cnt,Pfa);

% wall_x=save_det_data_dbscan(:,8);
% wall_y=save_det_data_dbscan(:,9);
% p = polyfit(wall_x,wall_y,1);
% x1 = linspace(-2,2);
% y1 = polyval(p,x1);
% figure
% plot(wall_x,wall_y,'o')
% hold on
% plot(x1,y1)
% hold off
% axis([-6, 6 0 10]);
% title('Wall Shape Regression')
% grid on

fclose(fid);