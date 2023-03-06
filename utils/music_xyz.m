function [x_vector, y_vector, z_vector] = music_xyz(Xcube,params,detout,num_tx,num_rx,azim_start,azim_end,azim_step,elev_start,elev_end,elev_step,n_target)
Nr=size(Xcube,1);   %%%length of Chirp
Ne=size(Xcube,2);   %%%length of channel: 12=4*3
Nd=size(Xcube,3);   %%%length of chirp loop

c = params.c; % Speed of light in air (m/s)
fc = params.fc; % Center frequency (Hz)
lambda = params.lambda;

num_detected_obj=size(detout,2);

virtual_ant=zeros(Ne,num_detected_obj); % 把每个接收天线存在物体的bin给取出来
for i=1:num_detected_obj
    virtual_ant(:,i)=Xcube(detout(2,i),:,detout(1,i));
end

angles=zeros(2,num_detected_obj); % each col: [azimuth, elevation]

% azimuth: -60~60, step: 1
azim_num=(azim_end-azim_start)/azim_step+1;
L = n_target;
a_azim = zeros(azim_num,2 * num_rx-L+1);
for al=0:1:(azim_num-1)
    for s_i = 1:1:(2 * num_rx-L+1)
        a_azim(al+1, s_i) = exp(1i*2*pi*fc*(s_i-1)*lambda/2*sin((azim_start+al*azim_step)/180*pi)/c);
    end
end
% elevation: -15~15, step: 1
elev_num=(elev_end-elev_start)/elev_step+1;
L = n_target;
a_elev = zeros(elev_num,2-L+1);
for al=0:1:(elev_num-1)
    for s_i = 1:1:(2-L+1)
        a_elev(al+1, s_i) = exp(1i*2*pi*fc*(s_i-1)*lambda/2*sin((elev_start+al*elev_step)/180*pi)/c);
    end
end

azimuth_ant = virtual_ant(1:2 * num_rx, :);
for i=1:num_detected_obj
    x=azimuth_ant(:,i);
    Rxx=x*x'; % TODO: 这里没有取平均吧？这样估计也太不准确了，多帧取平均，是不是会好些？现在是把每帧的点云图叠加
    e = eig(Rxx);
    [V, D] = eig(Rxx);
    Q = V(:,1:(size(Rxx,1)-n_target));
    test = Q*Q';
    maxv = -1;
    max_azim = 0;

    P_a = zeros(1, azim_num);
    for al=0:1:(azim_num-1)
        a_aa = a_azim(al+1, :);
        P_a(al+1) = 1/(a_aa * test * a_aa');
        if abs(P_a(al+1)) > maxv
            maxv = abs(P_a(al+1));
            max_azim = azim_start+al*azim_step;
        end
    end
    angles(1,i) =max_azim;
end


elevation_ant = virtual_ant([9,3], :);

for i=1:num_detected_obj
    x=elevation_ant(:,i);
    Rxx=x*x';
    e = eig(Rxx);
    [V, D] = eig(Rxx);
    Q = V(:,1:(size(Rxx,1)-n_target));
    test = Q*Q';
    maxv = -1;
    max_elev = 0;

    P_a = zeros(1, elev_num);
    for al=0:1:(elev_num-1)
        a_aa = a_elev(al+1, :);
        P_a(al+1) = 1/(a_aa * test * a_aa');
        if abs(P_a(al+1)) > maxv
            maxv = abs(P_a(al+1));
            max_elev = elev_start+al*elev_step;
        end
    end
    angles(2,i) =max_elev;
end

x_vector = cos(angles(2,:)*pi/180).*sin(angles(1,:)*pi/180);
y_vector = cos(angles(2,:)*pi/180).*cos(angles(1,:)*pi/180);
z_vector = sin(angles(2,:)*pi/180);

end