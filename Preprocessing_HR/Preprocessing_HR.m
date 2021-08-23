% The traget of this function is to preprocess the PPG signal using FIR
% band-pass filter and using Median filter to remove the baseline drift.
% Developer: Asher Wang
% Date: May 16, 2021

clear all
clc

load data_1.mat
time = round1.Times; % time
amplitude = round1.Channel2V;
Fs = 66.67;  % sampling frequency
figure;
plot(time,amplitude)

%% Kaiser Window Bandpass Filter
fcuts = [0.01 0.5 3 4];
mags = [0 1 0];
devs = [0.01 0.05 0.01];

[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,Fs);
n = n + rem(n,2);
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
[H,f] = freqz(hh,1,1024,Fs);

% plot the frequency response
figure;
plot(f,abs(H))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (abs)')

figure;
freqz(hh)  % plot in dB and rad/sample
grid on

%% FIR Bandpass Filter
sf1 = 0.3 * pi / Fs;
pf1 = 0.5 * pi / Fs;
pf2 = 3 * pi / Fs;
sf2 = 3.5 * pi / Fs;

pb = linspace(pf1,pf2,1e3)*pi;

%%%%%%%%%%FIR%%%%%%%%%%%%
bp = designfilt('bandpassfir', ...
    'StopbandAttenuation1',5, 'StopbandFrequency1',sf1,...
    'PassbandFrequency1',pf1,'PassbandRipple',0.01,'PassbandFrequency2',pf2, ...
    'StopbandFrequency2',sf2,'StopbandAttenuation2',30);

% bp = designfilt('bandpassfir','FilterOrder',40, ...
%          'CutoffFrequency1',0.5,'CutoffFrequency2',3, ...
%          'SampleRate',Fs);

%%%%%%%%%%%IIR%%%%%%%%%%%%
% bp = designfilt('bandpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',3, ...
%     'SampleRate',Fs);     

[h,w] = freqz(bp,1024);
hpb = freqz(bp,pb);

freqz(bp);

figure;
subplot(2,1,1)
plot(w/pi,abs(h),pb/pi,abs(hpb),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')

subplot(2,1,2)
plot(w/pi,db(h),pb/pi,db(hpb),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% Simple FIR Filter
b = fir1(20,[0.5*pi/Fs 3*pi/Fs]);
freqz(b,1,512)

%% Process the raw data using band-pass filter
signal_bandpass = filter(hh, 1, amplitude);
figure;
plot(time, signal_bandpass);
xlabel('Time (s)')
ylabel('Voltage (V)')
grid on

%% Median Filter to remove baseline drift
signal_median = medfilt1(signal_bandpass, 5);

figure;
plot(time, signal_bandpass, time, signal_median)
legend('Signal after band-pass filter', 'Signal after median filter')
xlabel('Time (s)')
ylabel('Voltage (V)')
grid on

%% Ellipord Filter (IIR zero phase shift digital filter to correct baseline drift)
wp = 3 * 2 / Fs;     % passband cutoff frequency
ws = 0.8 * 2 / Fs;     % stopband cutoff frequency 
pass_ripple = 0.005;    % passband ripple 
Rp = 20 * log10((1 + pass_ripple)/(1 - pass_ripple));   % passband ripple coefficient 
Rs = 20;                          % stopband attenuation 
[N Wn] = ellipord(Wp,Ws,Rp,Rs,'s');   % calculate the order of the elliptic filter 
[b a] = ellip(N,Rp,Rs,Wn,'high');       % calculate the coefficient of the elliptic filter 
[hw,w] = freqz(b,a,512);   
signal_ellipord = filter(b,a,signal_bandpass); 

figure
freqz(b,a);

figure
subplot(211); plot(time,signal_bandpass); 
xlabel('t(s)');ylabel('Amplitude');title('Signal after band-pass filter');grid
subplot(212); plot(time,result); 
xlabel('t(s)');ylabel('Amplitude');title('Signal after ellipord filter');grid
  
figure
N = 512
subplot(2,1,1);plot(abs(fft(signal_bandpass))*2/N);
xlabel('Frequency(Hz)');ylabel('Amplitude');title('Frequency response of the raw signal');grid;
subplot(2,1,2);plot(abs(fft(signal_ellipord))*2/N);
xlabel('Frequency(Hz)');ylabel('Amplitude');title('After linear filtering');grid;
ubplot(2,1,2);plot(abs(fft(signal_ellipord))*2/N);
xlabel('Signal spectrum after linear filtering');ylabel('Amplitude');grid;

%% HR Detect
figure
plot(time+60, result); grid on; % plot noisefreesignal 
[~,locs] = findpeaks(-signal_ellipord);  %get every peak point x-axis value, which back to locs
real_bottom_locs = intersect(find(signal_ellipord < 0), locs);  %only need value lager than 0.5, the max peak value
hold on
plot(real_bottom_locs/Fs, result(real_bottom_locs),'r*');   % add red star point at each peak in the figure
size = length(real_bottom_locs);                    % find the size of peak points (how many peak we have)
time1 = real_bottom_locs(size)/Fs - real_bottom_locs(1)/Fs;  % find the how many takes from first to last peak
HR = (size-1)/time1 * 60;     % first peak should counter 0 for the cycle and  calculate the heart rate

ylabel ('Voltage (V)')
xlabel ('Time (s)')
title ('Processed PPG signal')

str = {strcat('Heart Rate:' , num2str(HR), ' beats/min')} 
legend(str{:});   