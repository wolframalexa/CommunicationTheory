clear all; close all; clc

%% ECE300 Assignment 1
% Alexa Jakob
% September 14, 2020

% Question 1

% Set up filter
f = linspace(0,10,10000);
W = 1; % [dB]
H = exp(-log(2)/2 * (f/W).^2);
mag_H = 20*log10(abs(H)); % magnitude response

% d
T = 0.25;
A = 1/T;

% Fourier transform - spectrums
input_FT = A*T*exp(-j*pi*f*T).*sinc(f*T);
output_FT = input_FT.*H;
mag_output_FT = 20*log10(abs(output_FT));

% plot output
plot(f, mag_output_FT)
title('Output Magnitude Spectrum (Question 1d)')
xlabel('Frequency f (Hz)')
ylabel('Attenuation (dB)')
ylim([-60 0])


% e: use linear interpolation
linear_bounds = [-60 0];
B0 = interp1(mag_output_FT,f,-50); % reasonably linear between 2 sample points
attentuation_B0 = 20*log10(abs(A*T*exp(-j*pi*B0*T-log(2)/2*(B0/W).^2).*sinc(B0*T)));
% checked: for B0 = 3.3948, attenuation is -50dB and keeps descending

% g: graph output pulse
t = linspace(0, 4*T, 1000);
sigma = sqrt(log(2))/(2*pi*W);
y = A*(qfunc(-t/sigma)-qfunc((-t+T)/sigma));
figure
plot(t,y)
title('Output of rectangular pulse through Gaussian filter')
xlabel('Time')
ylabel('Amplitude')

%% Question 2
% Generate QAM grids
grid2 = GenerateQAMCoords(2);
grid4 = GenerateQAMCoords(4);
grid6 = GenerateQAMCoords(6);
grid8 = GenerateQAMCoords(8);

% compute epsilon_b on a decibel scale
k = 2;
eps_b2 = 1/k * 1/2.^k * sum(abs(grid2).^2);
eps_b2 = 10*log10(eps_b2);

k = 4;
eps_b4 = 1/k * 1/2.^k * sum(abs(grid4).^2);
eps_b4 = 10*log10(eps_b4);

k = 6;
eps_b6 = 1/k * 1/2.^k * sum(abs(grid6).^2);
eps_b6 = 10*log10(eps_b6);

k = 8;
eps_b8 = 1/k * 1/2.^k * sum(abs(grid8).^2);
eps_b8 = 10*log10(eps_b8);

% plot epsilons/dmin as a fcn of k
k = 2:2:8;
epsilons = [eps_b2 eps_b4 eps_b6 eps_b8];
figure
plot(k,epsilons)
xline(4)
xline(6)
title('\epsilon_b/d_{min}^2 vs k')
xlabel('k')
ylabel('\epsilon_b/d_{min}^2')

% calculate bits per symbol per dimension
etas = k/2;
figure
plot(k,etas, "*")
title('\eta: bits per symbol per dimension')
xlabel('k')
ylabel('\eta')
xlim([0 10])
ylim([0 5])
% Notice that the more symbols, the more bits per symbol per dimension.
% Although we usually want to minimize points in the constellation to
% reduce number of bits needed to be communicated, the larger the
% constellation, the greater the spectral efficiency.


% strategy for generating grids adapted from Nathaniel Kingsbury
function output = GenerateQAMCoords(k)
    n = 2.^(k/2);
    coords = -(n-1)/2:(n-1)/2;
    x_coords = repmat(coords, 1, n);
    y_coords = repelem(coords, n);
    output = complex(x_coords, y_coords);
end