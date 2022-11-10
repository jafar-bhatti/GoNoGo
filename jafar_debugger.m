%Make a pure tone
% fs = 192e3;
% dt = 1/fs; % seconds per sample 
% StopTime = 5; % seconds 
% t = (0:dt:StopTime)'; % seconds 
% F = 1000; % Sine wave frequency (hertz) 
% soundCard_output_scaling = 11; %Sound output from the sound card scaled up by 11. Divide by 11 to scale back down.
% data = sin(2*pi*F*t)/soundCard_output_scaling; 
% plot(t,data)
% ylim([-1.5, 1.5])
% data = [data, data]; %Make it have two channels
% PTB_data = data(:,1)';
%% Play sound using DAQ
fs = 192e3;
%Get stimulus
filt = load('booth1-221024-filter-192kHz.mat');
filt = filt.FILT;
[low, high] = generate_sound(filt);
figure 
spectrogram(high(:,1),256,200,256,fs,'yaxis')

% setup LYNX card (using DAQ)
d = daq.getDevices;
device = 'Lynx E44';
ch = [1 2];
description = sprintf('DirectSound Speakers (%s)',device);
ind = find(strcmp({d.Description},description));
dev = d(ind);
s = daq.createSession('directsound');
ch = addAudioOutputChannel(s,dev.ID,ch);
s.Rate = fs;
realFS = s.Rate;

%Queue output 
queueOutputData(s,high);
    
%Start output
startBackground(s);



%% Play and record pure tone using PsychToolBox
% playbackDevice = 'Speakers (Lynx E44)';
% recordingDevice = 'Record 01+02 (2- Lynx E44)'; % for booths 1-2, 3-4, device is 2-
% boothNumber = 1;        % Which booth we are calibrating, used to generate filter name
% 
% InitializePsychSound;
% devList = PsychPortAudio('GetDevices');
% windowsDSIdx = find(cell2mat(cellfun(@(X)~isempty(strfind(X,'MME')),{devList(:).HostAudioAPIName},'UniformOutput',false))); %#ok<STREMP>
% windowsDSIdx2 = find(cell2mat(cellfun(@(X)~isempty(strfind(X,'WASAPI')),{devList(:).HostAudioAPIName},'UniformOutput',false))); %#ok<STREMP>
% playbackIdx = find(cell2mat(cellfun(@(X)strcmp(X,playbackDevice),{devList(:).DeviceName},'UniformOutput',false)));
% recorderIdx = find(cell2mat(cellfun(@(X)strcmp(X,recordingDevice),{devList(:).DeviceName},'UniformOutput',false)));
% playbackIdx = intersect(playbackIdx,windowsDSIdx);
% recorderIdx = intersect(recorderIdx,windowsDSIdx);
% ph.player   = PsychPortAudio('Open',devList(playbackIdx).DeviceIndex,1,3,fs,1);
% ph.recorder = PsychPortAudio('Open',devList(recorderIdx).DeviceIndex,2,3,fs,1);
% 
% %record sound
% % Give the recorder an extra 2 seconds of buffer space over the player
% recorderBuffer = 2;
% PsychPortAudio('GetAudioData',ph.recorder,StopTime + recorderBuffer,StopTime + recorderBuffer);
% x.rec  = PsychPortAudio('Start',ph.recorder,1);
% 
% %play sound
% PsychPortAudio('FillBuffer',ph.player,PTB_data);
% x.play = PsychPortAudio('Start',ph.player,1);
% 
% %Wait for sound to end
% tic; WaitSecs(StopTime + recorderBuffer); toc
% 
% %Collect the recording
% [recNoise,~,~,x.recGet] = PsychPortAudio('GetAudioData',ph.recorder);
% PsychPortAudio('Stop',ph.recorder);
% 
% plot(recNoise*6) %Therefore, the sound input to the sound card is scaled down by 6. 
% 
% %Close PsychPort Audio
% PsychPortAudio('Close');

%% Play back my stimulus without filter and plot result 
playbackDevice = 'Speakers (Lynx E44)';
recordingDevice = 'Record 01+02 (2- Lynx E44)'; % for booths 1-2, 3-4, device is 2-
boothNumber = 1;        % Which booth we are calibrating, used to generate filter name
addpath(genpath('newFilters'));
addpath(genpath('_task'));
filt = load('booth1-221024-filter-192kHz.mat');
filt = filt.FILT;

%Constants
fs = 192e3;
rPa=20e-6;              % Refers to assumed air pressure in silence
vpPa=.316;              % Volts/Pascal conversion to get dB
inGain = 6;             % Sound card divides incoming stimulus by 6 
outGain = 11;           % Sound card multiplies outgoing stimulus by 11
lowerFreq = 3e3;        % Lower freq cutoff for filter
upperFreq = 70e3;       % Upper freq cutoff for filter (dB of filtered audio between low/upp should be ~equal)
noise_duration = .5;      % duration in seconds 

%MAKE STIMULUS

%Generate white noise and set spread
white_noise = randn(noise_duration*fs,1); %It's important that the array dimensions are the way they are for the convoluation and ramping
octave_spread = 1;
%Generate low frequency stim
mean1 = 10000;
[low1, high1] = octave_range(mean1, octave_spread);
[b, a] = butter(5, [low1, high1]/(fs/2)); %Generate a butterworth (bandpass) filter
low_stim = filtfilt(b, a, white_noise);
%Generate high frequency stim
mean2 = 30000;
[low2, high2] = octave_range(mean2, octave_spread);
[b, a] = butter(5, [low2, high2]/(fs/2)); %Generate a butterworth (bandpass) filter
high_stim = filtfilt(b, a, white_noise);


PTB_stim = (low_stim')/outGain;

%Start up PsychPort Audio
InitializePsychSound;
devList = PsychPortAudio('GetDevices');
windowsDSIdx = find(cell2mat(cellfun(@(X)~isempty(strfind(X,'MME')),{devList(:).HostAudioAPIName},'UniformOutput',false))); %#ok<STREMP>
windowsDSIdx2 = find(cell2mat(cellfun(@(X)~isempty(strfind(X,'WASAPI')),{devList(:).HostAudioAPIName},'UniformOutput',false))); %#ok<STREMP>
playbackIdx = find(cell2mat(cellfun(@(X)strcmp(X,playbackDevice),{devList(:).DeviceName},'UniformOutput',false)));
recorderIdx = find(cell2mat(cellfun(@(X)strcmp(X,recordingDevice),{devList(:).DeviceName},'UniformOutput',false)));
playbackIdx = intersect(playbackIdx,windowsDSIdx);
recorderIdx = intersect(recorderIdx,windowsDSIdx);
ph.player   = PsychPortAudio('Open',devList(playbackIdx).DeviceIndex,1,3,fs,1);
ph.recorder = PsychPortAudio('Open',devList(recorderIdx).DeviceIndex,2,3,fs,1);

% Record sound
recorderBuffer = .2;
PsychPortAudio('GetAudioData',ph.recorder,noise_duration + recorderBuffer,noise_duration + recorderBuffer); %preallocation
x.rec  = PsychPortAudio('Start',ph.recorder,1);

%play sound
PsychPortAudio('FillBuffer',ph.player,PTB_stim);
x.play = PsychPortAudio('Start',ph.player,1);

%Wait for sound to end and collect recorded sound 
tic; WaitSecs(noise_duration + recorderBuffer); toc
[recStim,~,~,t.recGet] = PsychPortAudio('GetAudioData',ph.recorder);
PsychPortAudio('Stop',ph.recorder);

%Plot unfiltered data 
noiseAdj = recStim(1,1000:end-1000) * inGain / rPa / vpPa;
[P,f] = pwelch(noiseAdj,1024,120,[],fs,'onesided');
dB = 10*log10(P);
f1 = figure(1); clf; 
hold on
plot(f,dB); drawnow;
disp(['Total volume ' num2str(10*log10(mean(P)*(f(end)-f(1))))...
    'dB in response to Jafar stim.']);

%% Apply filter and plot results
playbackDevice = 'Speakers (2- Lynx E44)';
recordingDevice = 'Record 01+02 (2- Lynx E44)'; % for booths 1-2, 3-4, device is 2-
boothNumber = 1;        % Which booth we are calibrating, used to generate filter name
addpath(genpath('newFilters'));
addpath(genpath('_task'));
filt = load('booth2-221024-filter-192kHz.mat');
filt = filt.FILT;

%Constants
fs = 192e3;
rPa=20e-6;              % Refers to assumed air pressure in silence
vpPa=.316;              % Volts/Pascal conversion to get dB
inGain = 6;             % Sound card divides incoming stimulus by 6 
outGain = 11;           % Sound card multiplies outgoing stimulus by 11
lowerFreq = 3e3;        % Lower freq cutoff for filter
upperFreq = 70e3;       % Upper freq cutoff for filter (dB of filtered audio between low/upp should be ~equal)
noise_duration = .5;      % duration in seconds 

%MAKE STIMULUS

%Generate white noise and set spread
white_noise = randn(noise_duration*fs,1); %It's important that the array dimensions are the way they are for the convoluation and ramping
octave_spread = 1;
%Generate low frequency stim
mean1 = 10000;
[low1, high1] = octave_range(mean1, octave_spread);
[b, a] = butter(5, [low1, high1]/(fs/2)); %Generate a butterworth (bandpass) filter
low_stim = filtfilt(b, a, white_noise);
%Generate high frequency stim
mean2 = 30000;
[low2, high2] = octave_range(mean2, octave_spread);
[b, a] = butter(5, [low2, high2]/(fs/2)); %Generate a butterworth (bandpass) filter
high_stim = filtfilt(b, a, white_noise);

%Add SAMs
low_stim_low_sam = add_SAM(low_stim,8,.9,3*pi/2,fs); %8Hz modulation
low_stim_high_sam = add_SAM(low_stim,16,.9,3*pi/2,fs); %16Hz modulation

high_stim_low_sam = add_SAM(high_stim,8,.9,3*pi/2,fs); %8Hz modulation
high_stim_high_sam = add_SAM(high_stim,16,.9,3*pi/2,fs); %16Hz modulation

%Filter noise
filtnoise_low = conv(low_stim_high_sam, filt, 'same'); %Convolve with speaker calibration filter
filtnoise_high = conv(high_stim_high_sam, filt, 'same'); %Convolve with speaker calibration filter
%ramp features
ramp_duration = 0.003;
%Add ramps
RT = ramp_duration*fs; %Ramp length in samples
ramp = (((-cos((0:(1/RT):1-1/RT)*pi)+1)/2).^2)'; %Create cosine squared ramp
ramper = ones(noise_duration*fs,1); %Create ramp enveloped
ramper(1:RT) = ramper(1:RT).*ramp; %Add ramp to beginning of ramp envelope
ramper(end-RT+1:end) = ramper(end-RT+1:end).*flipud(ramp); %Add ramp to end of ramp envelope
%Final stimuli
low_s = (filtnoise_low'.*ramper) ./outGain; %Apply ramp envelope to filtered noise and apply scaling factor
high_s = (filtnoise_high'.*ramper) ./outGain; %Apply ramp envelope to filtered noise and apply scaling factor

%flip
low_s = low_s';
high_s = high_s';

low_s = low_s .* 10^(10/20);

%where x is the number of decibels you want to
% increase or decrease (depending on if it is postive or negative).

high_s = high_s .* 10^(13/20);

%Start up PsychPort Audio
InitializePsychSound;
devList = PsychPortAudio('GetDevices');
windowsDSIdx = find(cell2mat(cellfun(@(X)~isempty(strfind(X,'MME')),{devList(:).HostAudioAPIName},'UniformOutput',false))); %#ok<STREMP>
windowsDSIdx2 = find(cell2mat(cellfun(@(X)~isempty(strfind(X,'WASAPI')),{devList(:).HostAudioAPIName},'UniformOutput',false))); %#ok<STREMP>
playbackIdx = find(cell2mat(cellfun(@(X)strcmp(X,playbackDevice),{devList(:).DeviceName},'UniformOutput',false)));
recorderIdx = find(cell2mat(cellfun(@(X)strcmp(X,recordingDevice),{devList(:).DeviceName},'UniformOutput',false)));
playbackIdx = intersect(playbackIdx,windowsDSIdx);
recorderIdx = intersect(recorderIdx,windowsDSIdx);
ph.player   = PsychPortAudio('Open',devList(playbackIdx).DeviceIndex,1,3,fs,1);
ph.recorder = PsychPortAudio('Open',devList(recorderIdx).DeviceIndex,2,3,fs,1);

% Record sound
recorderBuffer = .2;
PsychPortAudio('GetAudioData',ph.recorder,noise_duration + recorderBuffer,noise_duration + recorderBuffer); %preallocation
x.rec  = PsychPortAudio('Start',ph.recorder,1);

%play sound
PsychPortAudio('FillBuffer',ph.player,high_s);
x.play = PsychPortAudio('Start',ph.player,1);

%Wait for sound to end and collect recorded sound 
tic; WaitSecs(noise_duration + recorderBuffer); toc
[recStim,~,~,t.recGet] = PsychPortAudio('GetAudioData',ph.recorder);
PsychPortAudio('Stop',ph.recorder);

%Plot unfiltered data 
noiseAdj = recStim(1,1000:end-1000) * inGain / rPa / vpPa;
[P,f] = pwelch(noiseAdj,1024,120,[],fs,'onesided');

%look at the range that I care about. The other things in the power
%spectrum are from random things in the mic/room.
st = find(f>=5000,1,'first');
stop = find(f<=20000,1,'last');

dB = 10*log10(P);
f1 = figure(1); clf; 
hold on
plot(f,dB); drawnow;
title('Booth 2 Power Spectrum')
xlabel('Frequency')
ylabel('Decibels')
txt = ['Total volume ' num2str(10*log10(mean(P(st:stop))*(f(stop)-f(st)))) ' dB'];
text(5e4, 50, txt);
disp(txt)

%% Close audio devices
PsychPortAudio('Close'); %This is VERY VERY important
