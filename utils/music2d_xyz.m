function [x_vector, y_vector, z_vector] = music2d_xyz(Xcube,params,detout,num_tx,num_rx,azimuth_start,azimuth_end,azimuth_step,elevation_start,elevation_end,elevation_step,n_target)
Nr=size(Xcube,1);   %%%length of Chirp
Ne=size(Xcube,2);   %%%length of channel: 12=4*3
Nd=size(Xcube,3);   %%%length of chirp loop

c = params.c; % Speed of light in air (m/s)
fc = params.fc; % Center frequency (Hz)
lambda = params.lambda;
spacing=lambda/2;

num_detected_obj=size(detout,2);
range_indexes=detout(2,:);
doppler_indexes=detout(1,:);

virtual_ant=zeros(Ne,num_detected_obj);
for i=1:num_detected_obj
    virtual_ant(:,i)=Xcube(detout(2,i),:,detout(1,i));
end

angles=zeros(2,num_detected_obj); % each col: [azimuth, elevation]

% azimuth: -60~60, step: 1
azimuth_num=ceil((azimuth_end-azimuth_start)/azimuth_step)+1;
elevation_num=ceil((elevation_end-elevation_start)/elevation_step)+1;
L = n_target;
azimuth_sample=4;
elevation_sample=2;

a_theta_phi = zeros(azimuth_num,elevation_num,azimuth_sample,elevation_sample);
for theta=0:1:(azimuth_num-1)
    for phi=0:1:(elevation_num-1)
        for theta_i = 1:1:azimuth_sample
            for phi_i = 1:1:elevation_sample
                a_theta_phi(theta+1, phi+1,theta_i, phi_i) = exp(1i*2*pi*spacing/lambda*...
                ((theta_i-1)*sin(azimuth_start+theta*azimuth_step)*cos(elevation_start+phi*elevation_step)+...
            (phi_i-1)*cos(azimuth_start+theta*azimuth_step)*sin(elevation_start+phi*elevation_step)));
            end
        end
    end
end

for i=1:num_detected_obj
    x=virtual_ant(:,i);
    one_array=reshape(virtual_ant(:,i),4,3).'; % 3*4 [tx1; tx3; tx2]
    input_data = zeros(2,4);
    input_data(1,1) = one_array(3,1); % tx2 rx1
    input_data(1,2) = one_array(3,2); % tx2 rx2
    input_data(1,3) = one_array(3,3); % tx2 rx3
    input_data(1,4) = one_array(3,4); % tx2 rx4
    input_data(2,1) = one_array(1,3); % tx1 rx3
    input_data(2,2) = one_array(1,4); % tx1 rx4
    input_data(2,3) = one_array(2,1); % tx3 rx1
    input_data(2,4) = one_array(2,2); % tx3 rx2
    X=reshape(input_data.',[],1); %8*1
    Rxx = X*X';
    e = eig(Rxx);
    [V, D] = eig(Rxx);
    Q = V(:,1:(azimuth_sample*elevation_sample-n_target));
    test = Q*Q';
    maxv = -1;
    amaxd = 0;
    amaxr = 0;
    for theta=0:1:(azimuth_num-1)
        for phi=0:1:(elevation_num-1)
            aa_theta_phi = squeeze(a_theta_phi(theta+1, phi+1,:,:));
            aa_theta_phi=reshape(aa_theta_phi.',[],1).'; % 1*8
            P(theta+1, phi+1) = 1/(aa_theta_phi * test * aa_theta_phi');
            if abs(P(theta+1, phi+1)) > maxv
                maxv = abs(P(theta+1, phi+1));
                amaxd = theta;
                amaxr = phi;
            end
        end
    end
    angles(1,i) =azimuth_start+amaxd*azimuth_step;
    angles(2,i) =elevation_start+amaxr*elevation_step;
end


x_vector = cos(angles(2,:)*pi/180).*sin(angles(1,:)*pi/180);
y_vector = cos(angles(2,:)*pi/180).*cos(angles(1,:)*pi/180);
z_vector = sin(angles(2,:)*pi/180);

end