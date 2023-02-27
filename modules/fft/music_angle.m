function [AngData,agl_grid2] = music_angle(Xcube,params,ang_start,ang_end,ang_step,n_target)

Dopdata_mean = squeeze(mean(abs(Xcube), 2));
Nr=size(Xcube,1);   %%%length of Chirp
Ne=size(Xcube,2);   %%%length of receiver
Nd=size(Xcube,3);   %%%length of chirp loop

c = params.c; % Speed of light in air (m/s)
fc = params.fc; % Center frequency (Hz)
lambda = params.lambda;

% angle: -60~60, step: 1
ang_num=(ang_end-ang_start)/ang_step+1;
L = n_target;
a_a = zeros(ang_num,Ne-L+1);
for al=0:1:(ang_num-1)
    for s_i = 1:1:(Ne-L+1)
        a_a(al+1, s_i) = exp(1i*2*pi*fc*(s_i-1)*lambda/2*sin((ang_start+al*ang_step)/180*pi)/c);
    end
end


% win = taylorwin(Ne,5,-60);
% win = win/norm(win);
for i = 1:Nd
    for j = 1:Nr
        amp=Dopdata_mean(j,i);
        x=squeeze(Xcube(j,:,i));
        Rxx=x'*x;
        e = eig(Rxx);
        [V, D] = eig(Rxx);
        Q = V(:,1:(size(Rxx,1)-n_target));
        test = Q*Q';
        
        P_a = zeros(1, ang_num);
        for al=0:1:(ang_num-1)
            a_aa = a_a(al+1, :);
            P_a(al+1) = 1/(a_aa * test * a_aa');
        %         if abs(P(rl+1, theta+1)) > maxv
        %             maxv = abs(P(rl+1, theta+1));
        %             amaxd = rl;
        %             amaxr = theta;
        %         end
        end
        
       AngData(j,:,i)=P_a*amp;
    end
end

for i=0:1:(ang_num-1)
    agl_grid2(i+1)=ang_start+i*ang_step;
end
end