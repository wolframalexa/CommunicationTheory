clear; close all; clc

%% ECE300 Assignment 4
% Alexa Jakob
% November 2, 2020

%% Question 1: Graphing D(f||g)

sig1 = 1;
sig2 = linspace(0.1,10,1000);
Dfg = sig1^2./(2*sig2.^2) - 1/2;

plot(sig2,Dfg)
title("KL Divergence for continuous distributions f and g")
xlabel("\sigma_2")

yline(0, "--")
% Notice that the divergence dips below 0 for \sigma_2 ~ 1, which cannot
% happen in the discrete case

%% Question 3: Channel Capacity
p = linspace(0, 1/2, 1000);
% optimal prior
pi0 = 1./(1-3*p + (1-3*p).*exp((phi(p)-phi(2*p))./(1-3*p))) - 2*p./(1-3*p);
py0 = 2*p + pi0 - 3*p.*pi0; % P(Y = 0)
C = -pi0.*phi(p) - (1-pi0).*phi(2*p) + phi(py0); % channel capacity
C = 0.6932*C % convert to bits

% plot results
figure
plot(p,C)
hold on
plot(p,pi0)
title("Channel capacity and optimal prior")
legend("Channel capacity", "Optimal \pi_0", "location", "best")
xlabel("p")


%% Question 4: Water-Filling Algorithm

% Procedure:
% Pick a level \lambda
% if any component has \sigma_i^2 < \lambda, allocate no binits. distortion
% is \sigma_i^2
% else: allocate enough binits to reduce distortion from \sigma_i^2 to
% \lambda. distortion is lambda.

% to obtain R - D curve: R = R(\lambda), D = D(\lambda)

% parameters
sigmas = [1, 2, 3, 4];
n = 1000;
range = [0, 5];

% perform waterfilling algorithm
[tot_rate, tot_dist] = waterfill(sigmas, n, range);

% calculate entropy
H = sum(log(sqrt(2*pi*exp(1)*sigmas)));

% plot rate distortion curve
figure
plot(tot_dist, tot_rate)
title("Rate-Distortion Curve")
ylabel("Rate")
xlabel("Distortion")
hold on
plot(0, H, '*')
legend("R-D curve", "Vector Entropy")
% note that 0 distortion implies an infinte rate, although in theory the
% entropy should match R when D=0





%% Functions
function result = phi(p)
    % binary entropy function
    result = -p.*log(p)-(1-p).*log(1-p);
end


function [rate, dist] = waterfill(vars, numpts, lims)
% this function performs the water filling algorithm, given a vector of
% variances, a number of samples for lambda, and the limits of lambda

lambda = linspace(lims(1), lims(2), numpts);
rate = zeros(1,numpts);
dist = zeros(1,numpts);

numvars = size(vars,2);
tempdist = zeros(1,numvars);
for i = 1:numpts
    for j = 1:numvars
    % assign distortion to each component
        if lambda(i) <= vars(j)
            tempdist(j) = lambda(i);
        else
            tempdist(j) = vars(j);
        end
    end
    
    % find total distortion and rate
    dist(i) = sum(tempdist);
    rate(i) = sum(1/2 * log(vars./tempdist));
end
end

