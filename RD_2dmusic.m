%% compare range-doppler heatmap
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
Is_plot_rangeDop = 1;
Is_plot_rangeAng = 0;

data_path="D:\ti\mmwave_RawData\";
filename=data_path+"ego_motion\new_static3_Raw_0.bin";

%% 2D-FFT RANGE-DOPPLER HEATMAP
tic;
fid = fopen(filename,'r');
% data_frames=readDCA1000(filename);
for cnt=1:100
    disp(cnt)
%     data_frame = data_frames(:, (cnt-1)*data_each_frame+1:cnt*data_each_frame);
    data_frame=get_next_frame(fid,params);
    if cnt<15
        continue
    end
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
    Dopdata_sum = squeeze(abs(Dopplerdata_tx1(:,1,:)));
    
    % Plot range-Doppler image
    
    if Is_plot_rangeDop
        plot_rangeDop(Dopdata_sum,rng_grid,vel_grid);
    end
    break
end
fprintf(' \nIt took %6.2f s. \n',toc);
fclose(fid);

%% 2D-MUSIC RANGE-DOPPLER HEATMAP
tic;
fid = fopen(filename,'r');
for cnt=1:100
    data_frame=get_next_frame(fid,params);
    if cnt<15
        continue
    end
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
    break
end
toc

% =======================PARAMETERS=======================
n_target=10;
r_sample=120;
v_sample=60;
m_r=r_sample-n_target;
m_v=v_sample-n_target;
% range unit: m
range_start=0;
range_end=8;
range_step=0.02;
range_num=(range_end-range_start)/range_step+1;
% velocity unit: m/s
vel_start=-0.5;
vel_end=0.5;
vel_step=0.01;
vel_num=(vel_end-vel_start)/vel_step+1;

% =========================================================



chirp_tx1_avg = squeeze(chirp_tx1(:,1,:));
chirp_tx1_avg_crop=chirp_tx1_avg(1:r_sample,1:v_sample);

% range_sample_step=6;
% vel_sample_step=5;
% chirp_tx1_avg_crop=chirp_tx1_avg(1:range_sample_step:512,1:vel_sample_step:120);
% r_sample=size(chirp_tx1_avg_crop,1);
% v_sample=size(chirp_tx1_avg_crop,2);
% disp(r_sample)
% disp(v_sample)
m_r=r_sample-n_target;
m_v=v_sample-n_target;

% 2D spatial smoothing
l_r=r_sample-m_r+1;
l_v=v_sample-m_v+1;

X=zeros(m_r*m_v,l_r*l_v);
cnt2=0;
for p_r=1:l_r
    for p_v=1:l_v
        cnt2=cnt2+1;
        X_win=chirp_tx1_avg_crop(p_r:p_r+m_r-1,p_v:p_v+m_v-1);
        X_win=reshape(X_win,[],1);
        X(:,cnt2)=X_win;
    end
end

a_r = zeros(range_num,m_r);
% range: 0~10m, step: 1cm
for rl=0:1:(range_num-1)
    for s_i = 1:1:m_r
        a_r(rl+1, s_i) = exp(1i*4*pi*(range_start+rl*range_step)*sweepSlope/c*(s_i-1)*1/Fs);
%         a_r(rl+1, s_i) = exp(1i*4*pi*(range_start+rl*range_step)*sweepSlope/c*(1+range_sample_step*(s_i-1)-1)*1/Fs);
    end
end

% velocity: -2~2m/s, step: 0.01m/s
a_v = zeros(vel_num,m_v);
for vl=0:1:(vel_num-1)
    for s_i = 1:1:m_v
        a_v(vl+1, s_i) = exp(1i*4*pi*fc*(s_i-1)*Tc*(vel_start+vl*vel_step)/c);
%         a_v(vl+1, s_i) = exp(1i*4*pi*fc*(1+vel_sample_step*(s_i-1)-1)*Tc*(vel_start+vl*vel_step)/c);
    end
end
toc
% X=reshape(chirp_tx1_avg_crop,[],1);
J=exchange_matrix(m_r*m_v);
C = (X*X'+J*(X*X')'*J)/(2*l_r*l_v);
% e = eig(C);
% plot(e)
[V, D] = eig(C);
Q = V(:,1:(m_r*m_v-2*n_target));
test = Q*Q';
toc

P = zeros(range_num, vel_num);
maxv = -1;
amaxd = 0;
amaxr = 0;
for rl=0:1:(range_num-1)
    for vl=0:1:(vel_num-1)
        a_rr = a_r(rl+1, :);
        a_vv = a_v(vl+1, :);
        V_R_theta = kron(a_vv,a_rr);
        P(rl+1, vl+1) = 1/(V_R_theta * test * V_R_theta');
%         if abs(P(rl+1, theta+1)) > maxv
%             maxv = abs(P(rl+1, theta+1));
%             amaxd = rl;
%             amaxr = theta;
%         end
    end
end
toc
for i=0:1:(range_num-1)
    rng_grid2(i+1)=range_start+i*range_step;
end
for i=0:1:(vel_num-1)
    vel_grid2(i+1)=vel_start+i*vel_step;
end

P_crop=P(10:end,:);

figure(2)
[axh] = surf(vel_grid2,rng_grid2,abs(P));
view(0,90)
% set(gca,'YDir','normal')%对Y方向反转
colorbar
axis([-2 2 range_start range_end]);
grid off
shading interp

xlabel('Doppler Velocity (m/s)')
ylabel('Range(meters)')

% caxis([0,0.2])
title(['Range-Doppler heatmap by MUSIC (n\_target=',num2str(n_target),')'])

fprintf(' \nIt took %6.2f s. \n',toc);
fclose(fid);