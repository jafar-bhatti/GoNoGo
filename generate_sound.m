function [low, high] = generate_sound(filt)
%% set variables
fs = 192000;    % sample rate
noise_duration = 3;      % duration in seconds 

%Generate white noise and set spread
white_noise = randn(noise_duration*fs,1); %It's important that the array dimensions are the way they are for the convoluation and ramping
octave_spread = 1;

%Generate low frequency stim
mean1 = 10000;
[low1, high1] = range(mean1, octave_spread);
[b, a] = butter(5, [low1, high1]/(fs/2)); %Generate a butterworth (bandpass) filter
low_stim = filtfilt(b, a, white_noise);

%%Generate high frequency stim
mean2 = 30000;
[low2, high2] = range(mean2, octave_spread);
[b, a] = butter(5, [low2, high2]/(fs/2)); %Generate a butterworth (bandpass) filter
high_stim = filtfilt(b, a, white_noise);

%Filter noise
filtnoise_low = conv(low_stim, filt, 'same'); %Convolve with speaker calibration filter
filtnoise_high = conv(high_stim, filt, 'same'); %Convolve with speaker calibration filter

%ramp features
ramp_duration = 0.003;

%Add ramps
RT = ramp_duration*fs; %Ramp length in samples
ramp = (((-cos((0:(1/RT):1-1/RT)*pi)+1)/2).^2)'; %Create cosine squared ramp
ramper = ones(noise_duration*fs,1); %Create ramp enveloped
ramper(1:RT) = ramper(1:RT).*ramp; %Add ramp to beginning of ramp envelope
ramper(end-RT+1:end) = ramper(end-RT+1:end).*flipud(ramp); %Add ramp to end of ramp envelope

%Final stimuli
low_s = filtnoise_low.*ramper; %Apply ramp envelope to filtered noise
high_s = filtnoise_high.*ramper; %Apply ramp envelope to filtered nois

%Create chanhels channels
x = ones(length(low_s));
y = ones(length(high_s));

%Finals stimuli with channels
low = [low_s, x];
high = [high_s, y];






end 