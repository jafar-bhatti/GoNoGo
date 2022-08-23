function snd_mod = add_SAM(snd,mod_rate,depth,phase,fs)

% Add sinusoid amplitude modulation. 
% snd = sound to be modulated
% mod_rate = modulation rate
% depth = fraction of depth, 1=100%, 0.1=10%
% phase = starting phase of modulation (in radians) 3*pi/2 for zero
% fs = sample rate

% check sound scale
if round(mean(snd),1)~=0 % min(snd)>=0 - can't decide which if statement is better here
    snd = snd - mean(snd); % move the sound so the mean is at zero
end
% make the tone
n = 0:(length(snd))-1; % number of points
rps = 2*pi*mod_rate/fs; % radians per sample
x = n*rps + phase; % scale points by rps and add ph to change starting phase 
t = sin(x);
if depth>1
    depth = depth/100;
end
mod = (t/2)*depth+(1-depth/2);  % scale tone to between 0 and 1

if size(snd,1)>1
    snd = snd';
end
snd_mod = snd.*(mod);

