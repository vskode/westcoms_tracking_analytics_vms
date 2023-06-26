%% Frequency analysis
clear
close all;
% hold on
%% load
folder_path = '/home/vincent/Nicnak/third_paper/westcoms_trackingOutput/VMS/v09_3percent_vms/';
str_file = dir(folder_path);
last_file = load(join([folder_path, str_file(end).name]));
beached = last_file.tashore;
dbeached = diff(beached);

%% frequency parameters
sample_rate = 1/60;
L = numel(dbeached);

Y = fft(dbeached);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = sample_rate*(0:(L/2))/L;

%% find dominant frequency (tides)
[val, loc] = max(P1);
% the value of the wavelength corresponds to the tidal movement in hours
% its multiplied by 2 because it takes that long between two high tides
wavelength = 1/f(loc)/3600;

%% plot
semilogx(f,P1) 
title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f [Hz]")
ylabel("|P1(f)|")

%% remove high frequency influences from signal (i.e., wind, etc)

cutoff_frequency = 11;
Y_hp = fft(highpass(dbeached, 1/(cutoff_frequency*3600), 1/60));

T = ifft(fft(dbeached) - Y_hp);
figure()
plot(linspace(1, L, L)/3600, real(T))
title("Time resolved signal showing only tides and lower frequency components")
xlabel("t [hours]")
ylabel("number of newly beached particles")
