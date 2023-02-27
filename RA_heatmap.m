clc
clearvars
close all

addpath(genpath('./config'));
addpath(genpath('./utils'));
addpath(genpath('./modules/plot'));
addpath(genpath('./modules/fft'));
addpath(genpath('./modules/detection'));
addpath(genpath('./modules/cluster'));

% parameter setting
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
Is_plot_rangeDop = 0;
Is_plot_rangeAng = 1;

data_path="/home/agent/下载/";
% filename=data_path+"cornerRadar_1214\scenario4_right_Raw_0.bin";
filename=data_path+"radar1-0_person1_0_0_Raw_0.bin";

fid = fopen(filename,'r');

tic;
% data_frames=readDCA1000(filename);

M1=moviein(20);
filename1='test.gif';

for cnt=1:50
    disp(cnt)
%     data_frame = data_frames(:, (cnt-1)*data_each_frame+1:cnt*data_each_frame);
    data_frame=get_next_frame(fid,params);
%     if cnt<200
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
    Dopdata_sum = squeeze(mean(abs(Dopplerdata_tx1), 2));
    
%     if cnt==1
%         Dopdata_background=Dopdata_sum;
%     end
    % Plot range-Doppler image
    
    if Is_plot_rangeDop
        plot_rangeDop(Dopdata_sum,rng_grid,vel_grid);
    end
    
    Dopdata_merge=[Dopplerdata_tx1, Dopplerdata_tx3,Dopplerdata_tx2];
    
    
    % Angle FFT (Azimuth angle estimation)
    Angdata = fft_angle(Dopdata_merge(:,1:2 * Rx, :),fft_Ang,Is_Windowed);
    Angdata_dynamic=cat(3,Angdata(:,:,1:60),Angdata(:,:,62:120));
    
    [Angdata_music,agl_grid2] = music_angle(Dopdata_merge(:,1:2 * Rx, :),params,-60,60,1,1);
%     Angdata_music_dynamic=cat(3,Angdata_music(:,:,1:60),Angdata_music(:,:,62:120));
%     Angdata_music_static=Angdata_music(:,:,61);
    
    Angdata_music_dynamic=cat(3,Angdata_music(:,:,1:56),Angdata_music(:,:,66:120));
    Angdata_music_static=Angdata_music(:,:,57:65);
    
    if Is_plot_rangeAng
%         plot_rangeAng(Angdata_dynamic,rng_grid,agl_grid);
        plot_rangeAng(Angdata_music_dynamic,rng_grid,agl_grid2);
        plot_XY(Angdata_music_dynamic,rng_grid,agl_grid2,cnt);
    end
    
%     plot_rangeAng(Angdata_music,rng_grid,agl_grid2-30);
%     figure()
%     plot_XY(Angdata_music,rng_grid,agl_grid2-30,cnt);
    
    plot_XY(Angdata_music,rng_grid,agl_grid2-30,cnt);
    plot_rangeAng(Angdata_music,rng_grid,agl_grid2-30);
    M1(:,end+1)=getframe(gcf);
    [A,map]=rgb2ind(frame2im(getframe(gcf)),256);
    if cnt==1
        imwrite(A,map,filename1,'gif','Loopcount',inf,'DelayTime',0.3);
    else
        imwrite(A,map,filename1,'gif','WriteMode','append','DelayTime',0.3);
    end
    
    break
end

movie(M1,2)





