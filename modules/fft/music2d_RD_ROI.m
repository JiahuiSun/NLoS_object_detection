function [detout_music]=music2d_RD_ROI(params,chirp_avg,detout,n_target,r_sample,v_sample,delta_rng,delta_vel,range_step,vel_step)
detout_music=zeros(3,size(detout,2));% each column: [velocity, range, P]

m_r=r_sample-20;
m_v=v_sample-20;

% constant parameters
c = params.c; % Speed of light in air (m/s)
fc = params.fc; % Center frequency (Hz)

% configuration parameters
Fs = params.Fs;
sweepSlope = params.sweepSlope;

Tc = params.Tc; % us 

% Creat grid table
rng_grid = params.rng_grid;
vel_grid = params.vel_grid;
%% 2D spatial smoothing
chirp_avg_crop=chirp_avg(1:r_sample,1:v_sample);

l_r=r_sample-m_r+1;
l_v=v_sample-m_v+1;

X=zeros(m_r*m_v,l_r*l_v);
cnt2=0;
for p_r=1:l_r
    for p_v=1:l_v
        cnt2=cnt2+1;
        X_win=chirp_avg_crop(p_r:p_r+m_r-1,p_v:p_v+m_v-1);
        X_win=reshape(X_win,[],1);
        X(:,cnt2)=X_win;
    end
end
%% Correlation & SVD
J=exchange_matrix(m_r*m_v);
C = (X*X'+J*(X*X')'*J)/(2*l_r*l_v);
%     e = eig(C);
[V, D] = eig(C);
Q = V(:,1:(m_r*m_v-n_target));
test = Q*Q';

%% searching
for obj=1:size(detout,2)
    range_start=rng_grid(detout(2,obj))-delta_rng;
    range_end=rng_grid(detout(2,obj))+delta_rng;
    range_num=ceil((range_end-range_start)/range_step)+1;
    % velocity unit: m/s
    vel_start=vel_grid(detout(1,obj))-delta_vel;
    vel_end=vel_grid(detout(1,obj))+delta_vel;
    vel_num=ceil((vel_end-vel_start)/vel_step)+1;
    
    a_r = zeros(range_num,m_r);
    for rl=0:1:(range_num-1)
        for s_i = 1:1:m_r
            a_r(rl+1, s_i) = exp(1i*4*pi*(range_start+rl*range_step)*sweepSlope/c*(s_i-1)*1/Fs);
        end
    end

    a_v = zeros(vel_num,m_v);
    for vl=0:1:(vel_num-1)
        for s_i = 1:1:m_v
            a_v(vl+1, s_i) = exp(1i*4*pi*fc*(s_i-1)*Tc*(vel_start+vl*vel_step)/c);
        end
    end
    
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
            if abs(P(rl+1, vl+1)) > maxv
                maxv = abs(P(rl+1, vl+1));
                amaxd = rl;
                amaxr = vl;
            end
        end
    end  
    detout_music(2,obj)=range_start+amaxd*range_step;
    detout_music(1,obj)=vel_start+amaxr*vel_step;
    detout_music(3,obj)=maxv;
end

end