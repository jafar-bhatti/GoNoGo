function [low, high] = generate_sound_copy() %%% NOTE THIS WAS ADJUSTED TO OUTPUT UNFILTERED,UNRAMPED STIMULI
%% set variables
fs = 192000;    % sample rate
noise_duration = .5;      % duration in seconds 
event_duration = 0.01;    % duration in seconds
outGain = 11;             % output of sound card is multiplied by 11 so divide final stim by 11  

filt = load('practice_filter.mat');
filt = filt.FILT;

%Generate white noise and set spread
white_noise = randn(noise_duration*fs,1); %It's important that the array dimensions are the way they are for the convoluation and ramping
octave_spread = 1;

[P,f] = pwelch(white_noise,1024,120,[],fs,'onesided');
dB = 10*log10(P);
f1 = figure(1); clf; 
hold on
plot(f,dB); drawnow;
disp(['Total volume ' num2str(10*log10(mean(P)*(f(end)-f(1))))...
    'dB in response to flat noise.']);

figure
plot(white_noise)

%Generate low frequency stim
mean1 = 10000;
[low1, high1] = octave_range(mean1, octave_spread);
[b, a] = butter(5, [low1, high1]/(fs/2)); %Generate a butterworth (bandpass) filter
low_stim = filtfilt(b, a, white_noise);

[P,f] = pwelch(low_stim,1024,120,[],fs,'onesided');
dB = 10*log10(P);
f1 = figure(1); clf; 
hold on
plot(f,dB); drawnow;
title('bandPass_lowstim')
disp(['Total volume ' num2str(10*log10(mean(P)*(f(end)-f(1))))...
    'dB in response to flat noise.']);


%Generate high frequency stim
mean2 = 30000;
[low2, high2] = octave_range(mean2, octave_spread);
[b, a] = butter(5, [low2, high2]/(fs/2)); %Generate a butterworth (bandpass) filter
high_stim = filtfilt(b, a, white_noise);

%Add sinusoidal amplitude modulation (SAM)
low_stim_low_sam = add_SAM(low_stim,8,.9,3*pi/2,fs); %10Hz modulation
low_stim_high_sam = add_SAM(low_stim,16,.9,3*pi/2,fs); %20Hz modulation

high_stim_low_sam = add_SAM(high_stim,8,.9,3*pi/2,fs); %10Hz modulation
high_stim_high_sam = add_SAM(high_stim,16,.9,3*pi/2,fs); %20Hz modulation

%Plot sams
% figure
% subplot(2,2,1)
% spectrogram(low_stim_low_sam,256,200,256,fs,'yaxis')
% title('Low stim, Low amp mod') 
% 
% subplot(2,2,2)
% spectrogram(low_stim_high_sam,256,200,256,fs,'yaxis')
% title('Low stim, High amp mod') 
% 
% subplot(2,2,3)
% spectrogram(high_stim_low_sam,256,200,256,fs,'yaxis')
% title('High stim, Low amp mod') 
% 
% subplot(2,2,4)
% spectrogram(high_stim_high_sam,256,200,256,fs,'yaxis')
% title('High stim, High amp mod') 

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
low_s = (filtnoise_low.*ramper) ./outGain; %Apply ramp envelope to filtered noise and apply scaling factor
high_s = (filtnoise_high.*ramper) ./outGain; %Apply ramp envelope to filtered noise and apply scaling factor

%Apply decibel calculation (Kath)
%   stim = stim .* 10^(x/20) where x is the number of decibels you want to
%   increase or decrease (depending on if it is postive or negative).

%Length of event channel
eventChannelLength = fs*(noise_duration +(event_duration*2)); %add two events to the sound length

%Create event channel filled with zeros (no events)
x = zeros(eventChannelLength,1);
y = zeros(eventChannelLength,1);

%Add event at start and end to indicate sound onset and sound offset 
event = 0.5.*(ones(event_duration*fs,1));
x(1:length(event),1) = event;
x(((length(x)-length(event))+1):end, 1) = event;

y(1:length(event),1) = event;
y(((length(y)-length(event))+1):end,1) = event;

%Buffer sound with 0s so that the event channel and the sound channel are
%the same length (0.01s at each end)
pad = zeros(length(event), 1);
low_padded = [pad', low_s', pad'];
high_padded = [pad', high_s', pad'];

%Finals stimuli with channels
low = [low_padded', x];
high = [high_padded', y];

%Create time
% time = linspace(0,noise_duration+(2*event_duration),length(high));

%Plot the results (for testing only)
% figure
% subplot(3,1,1)
% spectrogram(low(:,1),256,200,256,fs,'yaxis')
% title('Low Frequency Bandwidth')
% subplot(3,1,2)
% plot(time, low(:,1))
% title('Sound Wave')
% subplot(3,1,3)
% plot(time, low(:,2))
% ylim([-.1 1])
% title('Event channel')
% 
% figure
% subplot(3,1,1)
% spectrogram(high(:,1),256,200,256,fs,'yaxis')
% title('High Frequency Bandwidth')
% subplot(3,1,2)
% plot(time, high(:,1))
% title('Sound Wave')
% subplot(3,1,3)
% plot(time, high(:,2))
% ylim([-.1 1])
% title('Event channel (0.5v signals event)')

end 