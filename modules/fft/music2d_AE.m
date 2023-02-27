function [angles_music]=music2d_AE(params,Dopdata_detout,azimuth_start,azimuth_end,azimuth_step,elevation_start,elevation_end,elevation_step)
lambda = params.lambda;
spacing=lambda/2;

angles_music=zeros(3,size(Dopdata_detout,1));% each column: [azimuth, elevation, P]
n_target=1;
azimuth_num=ceil((azimuth_end-azimuth_start)/azimuth_step)+1;
elevation_num=ceil((elevation_end-elevation_start)/elevation_step)+1;

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


for obj=1:size(Dopdata_detout,1)
    one_array=reshape(Dopdata_detout(obj,:),4,3).'; % 3*4 [tx1; tx3; tx2]
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
    C = X*X';
    %     e = eig(C);
    [V, D] = eig(C);
    Q = V(:,1:(azimuth_sample*elevation_sample-n_target));
    test = Q*Q';
    
    P = zeros(azimuth_num, elevation_num);
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
    angles_music(1,obj)=azimuth_start+amaxd*azimuth_step;
    angles_music(2,obj)=elevation_start+amaxr*elevation_step;
    angles_music(3,obj)=maxv;
    
end
end