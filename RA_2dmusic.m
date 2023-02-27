%% compare range-DOA heatmap
clc
clearvars
close all

addpath(genpath('.\config'));
addpath(genpath('.\utils'));
addpath(genpath('.\modules\plot'));
addpath(genpath('.\modules\fft'));
addpath(genpath('.\modules\detection'));
addpath(genpath('.\modules\cluster'));

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
Is_plot_rangeAng = 0;

data_path="D:\ti\mmwave_RawData\";
filename=data_path+"ego_motion\new_static3_Raw_0.bin";

%% 2D-FFT RANGE-ANGLE HEATMAP
fid = fopen(filename,'r');
% data_frames=readDCA1000(filename);
for cnt=1:100
    data_frame=get_next_frame(fid,params);
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
    
    % CFAR detector on Range-Velocity to detect targets 
    % Output format: [doppler index, range index(start from index 1), ...
    % cell power]
    Pfa = 1e-5; % probability of false alarm
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
            Rangedata_tx3(ri, :, :) = Rangedata_tx3(ri, :, :) * pha_comp_term_tx3;
            pha_comp_term_tx2 = exp(-1i * 4*pi * (detout(1,pick_idx)-fft_Vel/2-1) / (fft_Vel*3));
            Dopplerdata_tx2(ri, :, :) = Dopplerdata_tx2(ri, :, :) * pha_comp_term_tx2;
            Rangedata_tx2(ri, :, :) = Rangedata_tx2(ri, :, :) * pha_comp_term_tx2;
            
        end
    end
    
    Rangedata_merge = [Rangedata_tx1, Rangedata_tx3];
    Angdata = fft_angle(Rangedata_merge,fft_Ang,Is_Windowed);
    
    break
end
fclose(fid);
% Plot range-angle (RA) image
plot_rangeAng(Angdata,rng_grid,agl_grid);

%% 2D-MUSIC RANGE-ANGLE HEATMAP
a = zeros(121,8);
s = zeros(101,samples);

spacing=lambda/2;

% angle: 60~60, step: 1
for theta = 0:1:120
    for a_i=1:1:8
        a(theta+1, a_i) = exp(-1i*(a_i-1)*2*pi*cos((theta+30)/180*pi)*spacing/lambda);
    end
end

% range: 0~20m, step: 10cm
for rl=0:1:100
    for s_i = 1:1:samples
        s(rl+1, s_i) = exp(1i*4*pi*rl*sweepSlope/c*(s_i-1)/10*1/Fs);
    end
end

fid = fopen(filename,'r');
for cnt=1:100
    data_frame=get_next_frame(fid,params);
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
    
    chirp_tx1_avg = mean(chirp_tx1,3);
    chirp_tx3_avg = mean(chirp_tx3,3);
    chirp_merge_avg=[chirp_tx1_avg,chirp_tx3_avg];% 512*8
    
    X=reshape(chirp_merge_avg,[],1);
    n_target=1;
    
    C = X*X';
    e = eig(C);
    [V, D] = eig(C);
    Q = V(:,1:(512*8-n_target));
    test = Q*Q';
    
    P = zeros(101, 121);
    maxv = -1;
    amaxd = 0;
    amaxr = 0;
    for rl=0:1:100
        for theta = 0:1:120
            aa = a(theta+1, :);
            ss = s(rl+1, :);
            V_R_theta = kron(ss,aa);
            P(rl+1, theta+1) = 1/(V_R_theta * test * V_R_theta');
            if abs(P(rl+1, theta+1)) > maxv
                maxv = abs(P(rl+1, theta+1));
                amaxd = rl;
                amaxr = theta;
            end
        end
    end
    break
end
fclose(fid);
%% 
figure(2)
agl_grid2=zeros(1,121);
for i=1:1:121
    agl_grid2(i)=i-61;
end
rng_grid2=zeros(1,101);
for i=1:1:101
    rng_grid2(i)=i/10;
end
[axh] = surf(agl_grid2,rng_grid2,abs(P));
view(0,90)
% set(gca,'YDir','normal')%对Y方向反转
colorbar
axis([-60 60 0 10]);
grid off
shading interp

xlabel('Angle of arrive(degrees)')
ylabel('Range(meters)')

% caxis([0,0.2])
title('Range-Angle heatmap')

