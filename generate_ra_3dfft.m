%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used for processing the raw data collectd by TI awr1843
% radar
% Author : Xiangyu Gao (xygao@uw.edu), University of Washingyton
% Input: raw I-Q radar data
% Output: range-angle (RA) image, 3D point clouds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function  Angdata_all = generate_ra_3dfft(filename,domain)
clc;
clear;
close all;

addpath(genpath('.\config'));
addpath(genpath('.\utils'));
addpath(genpath('.\modules\plot'));
addpath(genpath('.\modules\fft'));
addpath(genpath('.\modules\detection'));
addpath(genpath('.\modules\cluster'));

% parameter setting
params = get_params_value();
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

% raw_data=readDCA1000(data_path+"sitting_posture_static_front_0223\"+action+".bin");

% Algorithm parameters
data_each_frame = samples*loop*Tx;
% set_frame_number = size(raw_data,2)/(data_each_frame)
% frame_start = 1;
% frame_end = set_frame_number;
Is_Windowed = 1;% 1==> Windowing before doing range and angle fft
Is_plot_rangeDop = 0;
Is_plot_rangeAng = 0;

Angdata_all=[];
point_cloud_feat_all=[];%frame num * 30*4

data_path="D:\ti\mmwave_RawData\";
action="center_d80_2";
filename=data_path+"sitting_posture_static_front_0223\"+action+".bin";
fid = fopen(filename);

% while ~feof(fid)
for cnt=1:100
%     disp(cnt)
    data_frame=get_next_frame(fid);
    % read the data of each frame, and then arrange for each chirps
%     data_frame = raw_data(:, (i-1)*data_each_frame+1:i*data_each_frame);
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
    Dopdata_sum = squeeze(mean(abs(Dopplerdata_tx1), 2));
    
    % Plot range-Doppler image
    if Is_plot_rangeDop
        plot_rangeDop(Dopdata_sum,rng_grid,vel_grid);
    end
    
    
    % CFAR detector on Range-Velocity to detect targets 
    % Output format: [doppler index, range index(start from index 1), ...
    % cell power]
    Pfa = 1e-4; % probability of false alarm
    [Resl_indx] = cfar_RV(Dopdata_sum, fft_Rang, num_crop, Pfa);
    detout = peakGrouping(Resl_indx);
    
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
    
    Rangedata_merge = [Rangedata_tx1, Rangedata_tx3];
    Dopdata_merge=[Dopplerdata_tx1, Dopplerdata_tx3,Dopplerdata_tx2];
    
    % Angle FFT (Azimuth angle estimation)
    Angdata = fft_angle(Dopdata_merge(:,1:2 * Rx, :),fft_Ang,Is_Windowed);
    [Angdata] = Normalize(Angdata, max_value);

    %% range-angle heatmap
    % Plot range-angle (RA) image
    if Is_plot_rangeAng
        plot_rangeAng(Angdata,rng_grid,agl_grid);
        saveas(gcf,'E:\大四下学期\毕设\figs_0223\side45\'+action+'_angle.jpg')
%         plot_rangeAng(Angdata_crop,rng_grid(num_crop+1:fft_Rang-num_crop),agl_grid);
    end
    
    
    Angdata_avg = abs(Angdata);
    Angdata_avg = squeeze(sum(Angdata_avg,3)/size(Angdata_avg,3)); %1024*128
    
    
    Angdata_avg_crop=Angdata_avg(1:40,:); %40*128 
    if domain==1
        Angdata_avg_trans=scaler(Angdata_avg_crop,rng_grid,agl_grid,0.6,0.8);%40*128 
    else
        Angdata_avg_trans=Angdata_avg_crop;
    end

    Angdata_avg_trans=reshape( Angdata_avg_trans,1,[]); %1*5120
    Angdata_all=[Angdata_all;Angdata_avg_trans];
    
    break
end
fclose(fid);






