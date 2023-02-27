function [x_vector, y_vector, z_vector] = naive_xyz(Xcube,detout,num_tx,num_rx,fft_Ang,Is_Windowed)
Nr=size(Xcube,1);   %%%length of Chirp
Ne=size(Xcube,2);   %%%length of channel: 12=4*3
Nd=size(Xcube,3);   %%%length of chirp loop

num_detected_obj=size(detout,2);
range_indexes=detout(2,:);
doppler_indexes=detout(1,:);

virtual_ant=zeros(Ne,num_detected_obj);
for i=1:num_detected_obj
    virtual_ant(:,i)=Xcube(detout(2,i),:,detout(1,i));
end


azimuth_ant = virtual_ant(1:2 * num_rx, :);
azimuth_fft=[];
for i=1:num_detected_obj
    if Is_Windowed
        win_xcube=azimuth_ant(:,i).*taylorwin(size(azimuth_ant,1));
    else
        win_xcube = azimuth_ant(:,i).*1;
    end
    azimuth_fft(:,i) =fft(win_xcube,fft_Ang);
end

[elem,k_max] = max(abs(azimuth_fft));  % shape = (1,num_detected_obj)
% peak_1 = azimuth_fft[k_max]
peak_1 = zeros(1,num_detected_obj);
for i =1:length(k_max)
    peak_1(i) = azimuth_fft(k_max(i), i);
end

for i=1:length(k_max)
    if k_max(i) > floor(fft_Ang / 2) - 1
        k_max(i)=k_max(i)-fft_Ang;
    end
end
wx = 2 * pi / fft_Ang * k_max;  % shape = (1,num_detected_obj)
x_vector = wx / pi;


elevation_ant = virtual_ant(2 * num_rx+1:Ne, :);
elevation_fft=[];
for i=1:num_detected_obj
    if Is_Windowed
        win_xcube=elevation_ant(:,i).*taylorwin(size(elevation_ant,1));
    else
        win_xcube = elevation_ant(:,i).*1;
    end
    elevation_fft(:,i) =fft(win_xcube,fft_Ang);
end

[elem,elevation_max] = max(abs(elevation_fft));
peak_2 = zeros(1,num_detected_obj);
for i =1:length(elevation_max)
    peak_2(i) = elevation_fft(elevation_max(i), i);
end

% Calculate elevation phase shift
wz = angle(peak_1 .* conj(peak_2) .* exp(1j * 2 * wx));
z_vector = wz / pi;
y_vector = sqrt(1 - x_vector .^ 2 - z_vector .^ 2);


end