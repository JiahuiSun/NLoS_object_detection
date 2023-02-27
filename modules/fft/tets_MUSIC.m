% clc;
% clear;
% close all;
addpath(genpath('E:\Project\mmWave\Workspace\point_cloud_demo\config'));
addpath(genpath('E:\Project\mmWave\Workspace\point_cloud_demo\utils'));
% addpath(genpath('.\modules\plot'));
% addpath(genpath('.\modules\fft'));
% addpath(genpath('.\modules\detection'));
% addpath(genpath('.\modules\cluster'));

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

% Algorithm parameters
data_each_frame = samples*loop*Tx;
% load pms1000_30fs.mat; 
data_path="D:\ti\mmwave_RawData\";
raw_data=readDCA1000(data_path+"human_position_0116\ver60_hori-60_2.bin");
% raw_data=readDCA1000(data_path+"wall_50.bin");
% read the data of each frame, and then arrange for each chirps
data_frame = raw_data(:, 1:data_each_frame);
data_chirp = [];
Tx=3;
loop=128;
for cj = 1:Tx*loop
    temp_data = data_frame(:, (cj-1)*samples+1:cj*samples);
    data_chirp(:,:,cj) = temp_data;
end
    
% separate the chirps for TDM-MIMO with 3 TXs
chirp_tx1 = data_chirp(:,:,1:3:end);
chirp_tx3 = data_chirp(:,:,2:3:end);
chirp_tx2 = data_chirp(:,:,3:3:end);
    
% permutation with the format [samples, Rx, chirp]
chirp_tx1 = permute(chirp_tx1, [2,1,3]);
chirp_tx2 = permute(chirp_tx2, [2,1,3]);
chirp_tx3 = permute(chirp_tx3, [2,1,3]);

chirp_merge=[chirp_tx1,chirp_tx3];

signal=chirp_tx1(:,:,1);
%%
NumElements=4;
c = physconst('LightSpeed');
fc = 77e9;              % Operating frequency
lambda = c/fc;
ElementSpacing=lambda/2;

Nsamp = 1024;
%%
%%
% Generate the multichannel signal received by the ULA.
%% Improving Resolution Using MVDR and MUSIC Estimators
NumSignals=1;
ScanAngles=[-90:90];


%% by WHJ 20180706 
[ymusic_my ang_my]= MUSIC(signal,NumSignals,NumElements,ScanAngles,1,ElementSpacing,fc,c);

figure(2)
PlotDOASpectra(ScanAngles,ymusic_my);
i=1;




function [scanpattern ang_estim]= MUSIC(X,NumSignals,NumElements,ang,L,elSp,freq,c)
%%% elsp: the interval between two antennas 
%%% NumElements: antenna number
%%% NumSignals: target number
%%% X: received signal
%%% L: subarray for spatial smoothing
%%% freq: center frequency
%%% c: speed
%%  From MUSICEstimator.m (privDOASpectrum function)
%% Compute eigenvectors of the covariance matrix
Nsig = NumSignals;   
NEle = NumElements;   
Cx=MLCovMtx(X,L);
% Cx=X.'*(X.')'/Nsig;
%% [eigenvals, eigenvects] = privEig(Cx);  
[eigenvects, eigenvalsC] = eig(Cx);
eigenvals = real(eigenvalsC);
[eigenvals,indx] = sort(diag(eigenvals),'descend');
eigenvects= eigenvects(:,indx);
eigenvals(eigenvals<0) = 0;
%% Form MUSIC denominator matrix from noise subspace eigenvectors
noise_eigenvects = eigenvects(:, Nsig + 1:end); 
%% position
EleIdx = 1:NEle;
delta = (NEle-1)/2+1;
numElIDX = numel(EleIdx);
pos  = [zeros(1,numElIDX);(EleIdx-delta)*elSp;zeros(1,numElIDX)];
%% tau returns the delay among sensor elements in a sensor array for a given direction specified in ANG. 
ang=[ang; zeros(1,length(ang))];
azang = ang(1,:);
elang = ang(2,:);
% angles defined in local coordinate system
incidentdir = [-cosd(elang).*cosd(azang);-cosd(elang).*sind(azang);-sind(elang)];
tau = pos.'*incidentdir/c;     

sv = exp(-1i*2*pi*freq*tau);    
D = sum(abs((sv'*noise_eigenvects)).^2,2)+eps(1); % 9.44 in Ref[1]
pPattern = 1./D; 
scanpattern = sqrt(pPattern);  
s = sign(diff(scanpattern));
iMax = 1+find(diff(s)<0);
[iPk iPk_pos]= sort(scanpattern(iMax),'descend');
iPk_pos=iMax(iPk_pos);
ang_estim=ang(1,iPk_pos(1:Nsig));
end

function Sx = MLCovMtx(X,L)
    % Maximum Likelihood
    K = size(X,1);
    M = size(X,2)-L+1;
    Sx = complex(zeros(M,M));
    for n=1:L
        % Forward-Only Spatial Smoothing
        subOut = X(:,n:n+M-1).';
        Sx =  Sx + 1/K*(subOut*subOut'); 
    end
    Sx = 1/L*Sx;
end
    
%% plot music
function PlotDOASpectra(x1,y1)
y1_dB = 20*log10(y1) - max(20*log10(y1));
plot(x1,y1_dB)
xlabel('Broadside Angle (degrees)');
ylabel('Power (dB)');
title('DOA Spatial Spectra')
legend('MUSIC');
grid on;

end





