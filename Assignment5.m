clear; close all; clc

%% ECE300 Assignment 5
% Alexa Jakob
% November 25, 2020

%% Question 1: Graphing Perr

gammas = linspace(5, 30, 10000);

mults = [20 12 8 8 4 4];
bk = [1 2 4 5 9 10];

% Calculate more accurate error probability
Perr_better = zeros(size(gammas));
for i=1:6
    element = 1/8 * mults(i) * qfunc(sqrt(bk(i)*gammas));
    Perr_better = Perr_better + element;
end

% Calculate less accurate error probability using min distance
Perr_est = 1/8 * 20 * qfunc(sqrt(gammas));

% plot both
gammas = 10*log10(gammas);
figure
semilogy(gammas, Perr_better);
hold on
semilogy(gammas, Perr_est);
title("BER Plot")
xlabel("\gamma_b, SNR/bit, dB")
ylabel("Probability of a bit error")
legend("Estimate with all terms", "Estimate with d_{min}", "location", "best")

%% Question 5: Simulating a Digital Communications System

% design filter
beta = 0.2; % rolloff factor
span = 4; % number of symbols
sps = 8; % samples per symbol
Rb = 1e6; % bits per second

a = rcosdesign(beta, span, sps);

% find impulse response
[p,t] = impz(a);

% g
q = conj(p(size(a,2)-t));
g = conv(q,p);
leng = size(g,1);
x = 1:leng; % support of g for graphing

% plot impulse response
figure
stem(t,p)
title("Impulse reponse")
xlabel("t")
ylabel("Amplitude")
xlim([0 size(t,1)])

% plot g
figure
stem(x,g)
title("g(t)")
xlabel("t")
ylabel("Amplitude")
xlim([0 leng])

% find peak of g
n0 = find(g==1);

% calculate worst case ISI
i = n0:8:leng;
i = i(2:end); % remove value @ n0
j = n0:-8:1;
j = j(2:end); % remove value @ n0
j = flip(j,2);
i = [j i]; % indices of ISI values
worstISI = sqrt(2)*sum(abs(g(i)));

sigpwr = 2; % QPSK: symbols all have mod sqrt(2)
SIR = 10*log10(sigpwr/worstISI^2);


% c: spectrum of g
W = (1+beta)* 1/2 * 1/2 * Rb;
Gf = fft(g,W);
magGf = abs(Gf);
f = 1/2 * Rb *(0:(W/2))/W;

% plot spectrum
figure
plot(f,magGf(1:W/2+1))
title("One sided spectrum of g")
xlabel("Frequency (Hz)")
ylabel("Amplitude")


%% d: simulation
% test with 1 + j
numbits = 1e6;
samplebits = [1, 1, zeros(1,numbits-2)];
[exBER, exmsdiff, exQPSK] = simulateQPSK(numbits, inf, samplebits, sps,p,q);

% test with various SNRs
reps = 10;
berinf = zeros(1,reps);
berminus5 = zeros(1,reps);
ber5 = zeros(1,reps);
ber10 = zeros(1,reps);

msinf = zeros(1,reps);
msminus5 = zeros(1,reps);
ms5 = zeros(1,reps);
ms10 = zeros(1,reps);

for i=1:reps
    bits = randi(2,[1 numbits])-1; % generate random independent bits
    
    [berinf(i), msinf(i)] = simulateQPSK(numbits,inf,bits,sps,p,q);
    [berminus5(i), msminus5(i)] = simulateQPSK(numbits,-5,bits,sps,p,q);
    [ber5(i), ms5(i)] = simulateQPSK(numbits,5,bits,sps,p,q);
    [ber10(i), ms10(i)] = simulateQPSK(numbits,10,bits,sps,p,q);    
end

% calculate means of BER
berinf_mean = sum(berinf)/reps;
berminus5_mean = sum(berminus5)/reps;
ber5_mean = sum(ber5)/reps;
ber10_mean = sum(ber10)/reps;
% BER values are in accordance with BER vs SNR curves given in the notes

% calculate means of mean squared diff
msinf_mean = sum(msinf)/reps;
msminus5_mean = sum(msminus5)/reps;
ms5_mean = sum(ms5)/reps;
ms10_mean = sum(ms10)/reps;

% calculate variances of BER
berinf_var = var(berinf);
berminus5_var = var(berminus5);
ber5_var = var(ber5);
ber10_var = var(ber10);

% calculate variances of mean squared diff
msinf_var = var(msinf);
msminus5_var = var(msminus5);
ms5_var = var(ms5);
ms10_var = var(ms10);

% collect results in a table
SNRvalues = [inf; -5; 5; 10];
MeanBER = [berinf_mean; berminus5_mean; ber5_mean; ber10_mean];
MeanMSDiff = [msinf_mean; msminus5_mean; ms5_mean; ms10_mean];
VarianceBER = [berinf_var; berminus5_var; ber5_var; ber10_var];
VarianceMSDiff = [msinf_var; msminus5_var; ms5_var; ms10_var];
T = table(SNRvalues, MeanBER, MeanMSDiff, VarianceBER, VarianceMSDiff)
% unsurprisingly, means & variances of BER and MSDiff decrease as SNR increases
% and the system becomes more precise

%% Function

function [BER, msdiff, QPSK] = simulateQPSK(numbits, gamma, bits, sps, p, q)
    % send signal
    odd = 1:2:size(bits,2);
    even = odd + 1;
    QPSK = 2*bits(odd)-1 + 2j * bits(even) - 1j; % generate QPSK
    
    % upsample 
    a_tilde = padarray(QPSK.',[0 sps-1],0,'post');
    a_tilde = reshape(a_tilde.', [1, numbits/2 * sps]);
    
    transmit = conv(p,a_tilde);
    
    % received symbol stream
    received = conv(transmit,q); % pass through matched filter
    
    % add noise
    noisepwr = 2*10^(-gamma/10);
    noise = sqrt(noisepwr)*(randn(1,size(received,2))+1j*randn(1,size(received,2))); % transform into standard normal
    received = received + noise;
    
    % downsample
    lenp = length(p);
    downsampleindex = lenp:8:(length(received)-lenp+1);
    sampled = received(1,downsampleindex);
    
    % decode
    decode_odd = double(real(sampled) > 0);
    decode_even = double(imag(sampled) > 0);
    decoded = reshape([decode_odd; decode_even], [1, numbits]);
    
    % compute BER
    BER = sum(abs(bits-decoded))/numbits;
    
    % compute mean square difference
    msdiff = sum(abs(QPSK-sampled).^2)/numbits;
end
