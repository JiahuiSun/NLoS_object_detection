function [elevAngdata] = elevation_angle_estim(Xcube,rng_grid,num_tx,num_rx,fft_Ang,Is_Windowed)
Nr=size(Xcube,1);   %%%length of Chirp
Ne=size(Xcube,2);   %%%length of channel: 12=4*3
Nd=size(Xcube,3);   %%%length of chirp loop

for i = 1:Nd
    for j = 1:Nr
        azimuth_ant = Xcube(j,1:2 * num_rx, i); %8*1
        if Is_Windowed
            win_xcube=azimuth_ant.*taylorwin(size(azimuth_ant,1));
        else
            win_xcube = azimuth_ant.*1;
        end
        azimuth_fft =fft(win_xcube,fft_Ang);
        [elem,k_max] = max(abs(azimuth_fft));
        peak_1 = azimuth_fft(k_max);
        if k_max > floor(fft_Ang / 2) - 1
            k_max=k_max-fft_Ang;
        end
        wx = 2 * pi / fft_Ang * k_max;
        x = wx / pi;
        
        elevation_ant = Xcube(j,2 * num_rx+1:Ne, i);%1*4
        if Is_Windowed
            win_xcube=elevation_ant.*taylorwin(size(elevation_ant,1));
        else
            win_xcube = elevation_ant.*1;
        end
        elevation_fft =fft(win_xcube,fft_Ang);
        [elem,elevation_max] = max(abs(elevation_fft));
        peak_2 = elevation_fft(elevation_max);


        % Calculate elevation phase shift
        wz = angle(peak_1 * conj(peak_2) * exp(1j * 2 * wx));
        z = wz / pi;
        R=rng_grid(j);
        if z/R>=-1 && z/R<=1
            elev_angle=asin(z/R)*180/pi;
        else
            elev_angle=-1000;
        end
        
        elevAngdata(j,i)=elev_angle;
    end
end
end