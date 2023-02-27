% plot range-angle heatmap
function [axh] = plot_angCurve(input,fft_Ang,agl_grid)
Ne = size(input,1);   %%%number of angleffts
c = physconst('LightSpeed');% Speed of light in air (m/s)
fc = 77e9; % Center frequency (Hz)

Is_Windowed=1;
if Is_Windowed
    win_xcube = input.*taylorwin(Ne);
else
    win_xcube = input;
end
AngData = fftshift(fft(win_xcube,fft_Ang));

ScanAngles=[-90:90];
[scanpattern,ang_estim]= MUSIC(input,1,8,ScanAngles,1,c/(fc*2),fc,c);
figure('visible','on')
subpplot(121)
[axh] = plot(agl_grid,abs(AngData))
xlim([-90 90]);
% grid off
% shading interp
xlabel('Angle of arrive(degrees)')
ylabel('amp')

subpplot(122)
[axh] = plot(ang_estim)

title('Azimuth-Angle estimation')
end