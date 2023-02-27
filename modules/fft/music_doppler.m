function [DopData]=music_doppler(Xcube,params,vel_start,vel_end,vel_step,n_target)
%%%%  Xcube : Nr*Ne*Nd , original data
%%%   fft_Rang: range fft length
%%%   fft_Vel:  velocity fft length(2D-FFT)
%%%   fft_Ang:  angle fft length(3D FFT)

Nr=size(Xcube,1);   %%% length of Chirp
Ne=size(Xcube,2);   %%% # of receiver
Nd=size(Xcube,3);   %%% # of chirp loop

c = params.c; % Speed of light in air (m/s)
fc = params.fc; % Center frequency (Hz)
Tc = params.Tc; % us

% range_amp=abs(Xcube(:,1,1));

%% Second fft on dopper dimension
% velocity: -2~2m/s, step: 0.01m/s

vel_num=(vel_end-vel_start)/vel_step+1;
L = n_target;
a_v = zeros(vel_num,Nd-L+1);
for vl=0:1:(vel_num-1)
    for s_i = 1:1:(Nd-L+1)
        a_v(vl+1, s_i) = exp(1i*4*pi*fc*(s_i-1)*Tc*(-1)*(vel_start+vl*vel_step)/c);
    end
end
for i=1:Ne
    for j=1:Nr
        x=squeeze(Xcube(j,i,:));
        
        Rxx=x*x';
        RSM = spsmooth(Rxx,L);
        e = eig(RSM);
        [V, D] = eig(RSM);
        Q = V(:,1:(size(RSM,1)-n_target));
        test = Q*Q';

        P_v = zeros(1, vel_num);
        for vl=0:1:(vel_num-1)
            a_vv = a_v(vl+1, :);
            P_v(vl+1) = 1/(a_vv * test * a_vv');
        %         if abs(P(rl+1, theta+1)) > maxv
        %             maxv = abs(P(rl+1, theta+1));
        %             amaxd = rl;
        %             amaxr = theta;
        %         end
        end
        
       DopData(j,i,:)=P_v;
     end
end
end