%% ECE300 Assignment 3
% Alexa Jakob
% October 5, 2020

clear; close all; clc

%% 2
% parameters of question
W = 10000; % message bandwidth, Hertz
D = [0.2 1 5]; % dispersion coefficients
A = 2; % peak amplitude

Df = D*W; % frequency derivation
kf = Df/A; % frequency sensitivity
bw_carson = 2*Df + 2*W*ones(1,3); % bandwidth estimated with carson's rule
% betas ~= D
bw_curve = [15 6 3].*Df; % by inspection of curve

%% 3

%% AM
alpha = 0.3;
t0 = 5;
T = 10;
fs = 10000; % sampling rate
N = fs*T + 1;
t = linspace(0,T,N);

m = alpha./(alpha^2+(t-t0).^2);
mhat = m/alpha.*(t-t0);

m0 = m - mean(m); % subtract mean for easier math
maxm0 = max(abs(m0));

% a
fc = 20; % Hz

% baseband signal is m0(t)

% bandpass signals
% envelope = sqrt(i^2 + q^2)
% AM with modulation index 80%
ka = maxm0/0.8;
i_AM = 1+ka*m0;
x_AM = i_AM.*cos(2*pi*fc*t); % assume carrier amplitude 1
env_AM = sqrt(i_AM.^2);

% DSB-SC
x_DSBSC = m0.*cos(2*pi*fc*t);
env_DSBSC = sqrt(m0.^2);

% USSB
x_USSB = m0.*cos(2*pi*fc*t) - mhat.*sin(2*pi*fc*t);
env_USSB = sqrt(m0.^2 + mhat.^2);

% LSSB
i_LSSB = m0;
q_LSSB = -mhat;
x_LSSB = i_LSSB.*cos(2*pi*fc*t) - q_LSSB.*sin(2*pi*fc*t);
env_LSSB = sqrt(i_LSSB.^2 + q_LSSB.^2);

% plot message and hilbert transform
plot(t,m0)
hold on
plot(t,mhat)
xlabel("Time")
ylabel("Message Amplitude")
legend("Message", "Hilbert Trans.")
title("Message and Hilbert Transform")
hold off

% calculate in dB
N0 = 2^17;
M0_f = fft(m0, N0);
abs_M0f = 20*log10(abs(M0_f));
abs_M0f = abs_M0f(1:N0/2+1);

Mhat_f = fft(mhat, N0);
abs_Mhatf = 20*log10(abs(Mhat_f));
abs_Mhatf = abs_Mhatf(1:N0/2+1);

f = fs*(0:(N0/2))/N0;

% plot dB calculations
figure
plot(f,abs_M0f)
hold on
plot(f,abs_Mhatf)
title("Response in dB")
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)")
legend("Message", "Hilbert trans.","Location","best")
xlim([0 10])

% plot for each format
figure
subplot(2,2,1)
plot(t,x_AM)
hold on
plot(t,env_AM)
title("AM")
xlabel("Time")

subplot(2,2,2)
plot(t,x_DSBSC)
hold on
plot(t,env_DSBSC)
title("DSB-SC")
xlabel("Time")

subplot(2,2,3)
plot(t,x_USSB)
hold on
plot(t,env_USSB)
title("USSB")
xlabel("Time")

subplot(2,2,4)
plot(t,x_LSSB)
hold on
plot(t,env_LSSB)
title("LSSB")
xlabel("Time")

% calculate magnitude spectra for each format
FT_AM = fft(x_AM,N0);
abs_FT_AM = 20*log10(abs(FT_AM));
abs_FT_AM = abs_FT_AM(1:N0/2+1);

FT_DSBSC = fft(x_DSBSC,N0);
abs_FT_DSBSC = 20*log10(abs(FT_DSBSC));
abs_FT_DSBSC = abs_FT_DSBSC(1:N0/2+1);

FT_USSB = fft(x_USSB,N0);
abs_FT_USSB = 20*log10(abs(FT_USSB));
abs_FT_USSB = abs_FT_USSB(1:N0/2+1);

FT_LSSB = fft(x_LSSB,N0);
abs_FT_LSSB = 20*log10(abs(FT_LSSB));
abs_FT_LSSB = abs_FT_LSSB(1:N0/2+1);

% plot magnitude spectra
figure
subplot(2,2,1)
plot(f,abs_FT_AM)
title("AM")
xlabel("Frequency (Hz)")
ylabel("Spectrum (dB)")
xlim([0 50])

subplot(2,2,2)
plot(f,abs_FT_DSBSC)
title("DSB-SC")
xlabel("Frequency (Hz)")
ylabel("Spectrum (dB)")
xlim([0 50])

subplot(2,2,3)
plot(f,abs_FT_USSB)
title("USSB")
xlabel("Frequency (Hz)")
ylabel("Spectrum (dB)")
xlim([0 50])

subplot(2,2,4)
plot(f,abs_FT_LSSB)
title("LSSB")
xlabel("Frequency (Hz)")
ylabel("Spectrum (dB)")
xlim([0 50])

%% FM
kf1 = 1;
kf2 = 10;
Ts = 1/fs;
intm0 = cumsum(Ts*m0);

% Time Domain (fc = 20Hz)
fc = 20;
x_FM_fc1_kf1 = cos(2*pi*fc*t+2*pi*kf1*intm0);
x_FM_fc1_kf2 = cos(2*pi*fc*t+2*pi*kf2*intm0);

% plot FM signals
figure
subplot(2,1,1)
plot(t,x_FM_fc1_kf1)
xlim([4 6])
title('kf = 1')
xlabel("Time (s)")

subplot(2,1,2)
plot(t,x_FM_fc1_kf2)
xlim([4 6])
title('kf = 10')
xlabel("Time (s)")

sgtitle('Time Domain (fc = 20 Hz)')

% Frequency Domain (fc = 1kHz)
fc = 1000;
x_FM_fc2_kf1 = cos(2*pi*fc*t + 2*pi*kf1*intm0);
x_FM_fc2_kf2 = cos(2*pi*fc*t + 2*pi*kf2*intm0);

FT_FM1 = fft(x_FM_fc2_kf1,N0);
abs_FT_FM1 = 20*log10(abs(FT_FM1));
abs_FT_FM1 = abs_FT_FM1(1:N0/2+1);

FT_FM2 = fft(x_FM_fc2_kf2,N0);
abs_FT_FM2 = 20*log10(abs(FT_FM2));
abs_FT_FM2 = abs_FT_FM2(1:N0/2+1);

figure
subplot(2,1,1)
plot(f,abs_FT_FM1)
axis([fc-10, fc+80, 0, 90])
title('kf = 1')
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)")

subplot(2,1,2)
plot(f,abs_FT_FM2)
axis([fc-10, fc+80, 0, 90])
title('kf = 10')
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)")

sgtitle("Magnitude Spectrum of FM signals")


% bandwidth analysis
W = 3; % Hz

df1 = kf1*maxm0;
df2 = kf2*maxm0; % frequency deviations

% use Carson's rule
BW1 = 2*df1 + 2*W; % 12.0623
BW2 = 2*df2 + 2*W; % 66.6233

% This makes sense - increasing kf tenfold yields an increase 10/2 =
% 5-fold. It is also reflected on the graphs' bandwidth.

% The FM spectrum is not comprised of spectral lines since m0(t) is not a
% pure tone; however, we can see "bounces" in the spectrum for the length
% of the bandwidth - there are still peaks in the magnitude spectrum, and
% increasing kf increases the number and width of the peaks.