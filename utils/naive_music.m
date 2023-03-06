function [x_vector, y_vector, z_vector,RCS] = naive_music(Xcube,detout,num_tx,num_rx,fft_Ang,Is_Windowed)
Nr=size(Xcube,1);   %%%length of Chirp
Ne=size(Xcube,2);   %%%length of channel: 12=4*3
Nd=size(Xcube,3);   %%%length of chirp loop

num_detected_obj=size(detout,2);

virtual_ant=zeros(Ne,num_detected_obj);
RCS=zeros(num_detected_obj);
for i=1:num_detected_obj
    virtual_ant(:,i)=Xcube(detout(2,i),:,detout(1,i));
    RCS(i)=sum(abs(Xcube(detout(2,i),:,detout(1,i))));
end

RCS=(RCS-min(RCS))/max(RCS);

azimuth_ant = virtual_ant(1:2 * num_rx, :);
azimuth_music=[];
for i=1:num_detected_obj
    [values,~]=pmusic(azimuth_ant(:,i),1);
    azimuth_music(:,i) = values;
end

[~,k_max] = max(abs(azimuth_music));  % shape = (1,num_detected_obj)
disp(["Here" length(azimuth_music(1,:)) length(azimuth_music(:,1)) length(k_max)]);
for i=1:length(k_max)
    k_max(i)=k_max(i)/256*2;
    if k_max(i)>1
        k_max(i)=k_max(i)-2;
    end
end
x_vector = k_max;

ks=zeros([4 num_detected_obj]);
elevation_music=zeros([256 num_detected_obj]);
for j=4
    elevation_ant = virtual_ant([j+8 j+2], :);
    for i=1:num_detected_obj
        [values,~]=pmusic(elevation_ant(:,i),1);
        elevation_music(:,i) = values;
    end
    [~,k2_max] = max(abs(elevation_music));  % shape = (1,num_detected_obj)
    ks(j,:)=k2_max;
end
% figure(100);
% scatter(real(elevation_ant(1,:)),imag(elevation_ant(1,:)));
% hold on
% scatter(real(elevation_ant(2,:)),imag(elevation_ant(2,:)));
% hold on
% for i=1:num_detected_obj
%     plot([real(elevation_ant(1,i)),real(elevation_ant(2,i))],[imag(elevation_ant(1,i)),imag(elevation_ant(2,i))],Color='b')
% save('ks.mat','ks');
% ks=mean(ks);
% disp(["ks" ks]);
% disp(["k2" k2_max]);
% for i=1:length(ks)
%     ks(i)=ks(i)/256*2;
%     if ks(i)>1
%         ks(i)=ks(i)-2;
%     end
% end
% z_vector=ks;

for i=1:length(k2_max)
    k2_max(i)=k2_max(i)/256*2;
    if k2_max(i)>1
        k2_max(i)=k2_max(i)-2;
    end
end

z_vector=k2_max;

% Calculate elevation phase shift

y_vector = sqrt(1 - x_vector .^ 2 - z_vector .^ 2);
disp(["music" x_vector]);
end