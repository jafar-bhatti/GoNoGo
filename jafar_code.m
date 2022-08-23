%% set variables
fs = 200000;    % sample rate
dur = 0.3;      % duration in seconds (300 ms)

%% Determine filter range

%This code generates different octave spreads. The first value in the low
%array is associated with the first value is the high array. For example, 
%the difference between low(1) and high(1) is 1-octave spread around the 
%mean frequency 18,000 kHz. The difference between low(9) and high(9) is a
%2-octave spread around the mean frequency 18,00 kHz. 

mean_freq = 18000; %The guassian is centered around this frequency
octaves = 1:(0.25/2):2; %Add 1 to 2 octaves in each direction, iterate by 0.125 octaves. Adding 1 octave in each direction means the spread is 2 octaves.
[low, high] = range(18000, octaves); %Extract the max and min frequency values for each octave shift 

%% Bandwidth Modulation
white_noise = randn(1,fs); %Generate white noise (Seed selection???)

%2-octave spread 
[b, a] = butter(5, [low(1), high(1)]/(fs/2)); %Generate a butterworth (bandpass) filter
low_base = filtfilt(b, a, white_noise); %Apply filter
spectrogram(low_base,256,200,256,fs,'yaxis'); %View the stimulus
title("Gaussian White Noise (2 Octave spread)");

%4-octave spread
[b, a] = butter(5, [low(9), high(9)]/(fs/2)); %Generate a butterworth (bandpass) filter
high_base = filtfilt(b, a, white_noise); %Apply filter
figure
spectrogram(high_base,256,200,256,fs,'yaxis'); %View the stimulus
title("Gaussian White Noise (4 Octave spread)");

%Gaussian noise with increasing octave spread 
figure
tiledlayout(2, 4);
for i = 1:8
[b, a] = butter(5, [low(i+1), high(i+1)]/(fs/2));
y = filtfilt(b, a, white_noise); %Apply filter
nexttile;
spectrogram(y,256,200,256,fs,'yaxis'); %View the stimulus
title("Bandwidth = " + octaves(i+1)*2 + " Octaves");
end

%% Amplitude Modulation
mod_rate = [3, 5, 10];
figure
for i = 1:length(mod_rate)
mod_noise = add_SAM(low_base,mod_rate(i),.9,3*pi/2,fs);
subplot(3, 1, i);
spectrogram(mod_noise,256,200,256,fs,'yaxis');
title("Amp Mod Rate = " + mod_rate(i));
end

%% Mean shifted Gaussian white noise
white_noise = randn(1,fs);
octave_spread = 1;

mean1 = 10000;
[low1, high1] = range(mean1, octave_spread);
[b, a] = butter(5, [low1, high1]/(fs/2)); %Generate a butterworth (bandpass) filter
low_stim = filtfilt(b, a, white_noise);
figure
spectrogram(low_stim,256,200,256,fs,'yaxis')
title("Low frequency Gaussian white noise (μ = 10kHz)");

mean2 = 30000;
[low2, high2] = range(mean2, octave_spread);
[b, a] = butter(5, [low2, high2]/(fs/2)); %Generate a butterworth (bandpass) filter
high_stim= filtfilt(b, a, white_noise);
figure
spectrogram(high_stim,256,200,256,fs,'yaxis')
title("High frequency Gaussian white noise (μ = 30kHz)");

%% Ampltiude shifted Gaussian white noise
low_amp = 3;
high_amp = 10;

low_stim_low_amp = add_SAM(low_stim,low_amp,.9,3*pi/2,fs);
figure
spectrogram(low_stim_low_amp,256,200,256,fs,'yaxis')
title("Low frequency white noise w/ low amp mod");

low_stim_high_amp = add_SAM(low_stim,high_amp,.9,3*pi/2,fs);
figure
spectrogram(low_stim_high_amp,256,200,256,fs,'yaxis')
title("Low frequency white noise w/ high amp mod");

high_stim_low_amp = add_SAM(high_stim,low_amp,.9,3*pi/2,fs);
figure
spectrogram(high_stim_low_amp,256,200,256,fs,'yaxis')
title("High frequency white noise w/ low amp mod");

high_stim_high_amp = add_SAM(high_stim,high_amp,.9,3*pi/2,fs);
figure
spectrogram(high_stim_high_amp,256,200,256,fs,'yaxis')
title("High frequency white noise w/ high amp mod");

%% Functions
function [low, high] = range(mean_freq, shift)

low = 2.^(log2(mean_freq) - shift); 
high = 2.^(log2(mean_freq) + shift);

end