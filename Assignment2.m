%% ECE300 Assignment 2
% Alexa Jakob
% September 25, 2020

clear all; close all; clc

% Question 1
% a
r1 = 1; % number of red balls in urn 1
b1 = 9; % number of blue balls in urn 1
r2 = 1; % number of red balls in urn 2
b2 = 9; % number of blue balls in urn 2

% prior distribution - probability of getting red/blue from urn 1
pi_R = r1/(r1+b1);
pi_B = b1/(r1+b1);

% likelihood fcn - given result from urn 1, prob of pulling from urn 2
Plikely_R_B = r2/(r2+b2+1);
Plikely_B_B = (b2+1)/(r2+b2+1);
Plikely_R_R = (r2+1)/(r2+b2+1);
Plikely_B_R = b2/(r2+b2+1);

bayes_sum_R = Plikely_R_B*pi_B + Plikely_R_R*pi_R;
bayes_sum_B = Plikely_B_B*pi_B + Plikely_B_R*pi_R;

% posteriori distribution - given urn2, return urn1
Ppost_R_B = Plikely_B_R*pi_R/bayes_sum_B;
Ppost_R_R = Plikely_R_R*pi_R/bayes_sum_R;
Ppost_B_R = Plikely_R_B*pi_B/bayes_sum_R;
Ppost_B_B = Plikely_B_B*pi_B/bayes_sum_B;

% b & c
% make decision vectors
MAP = [0 0];
ML = [0 0];

% MAP maximizes a-posteriori distribution
if Ppost_B_R > Ppost_R_R
    MAP(1) = 2;
    MAP_error = bayes_sum_R * Ppost_R_R;
else
    MAP(1) = 1;
    MAP_error = bayes_sum_R * Ppost_B_R;
end

if Ppost_B_B > Ppost_R_B
    MAP(2) = 2;
    MAP_error = MAP_error + bayes_sum_B * Ppost_R_B;
else
    MAP(2) = 1;
    MAP_error = MAP_error + bayes_sum_B * Ppost_B_B;
end

%MAP_error = 1 - MAP_error;

% ML maximizes likelihood fcn
if Plikely_R_B > Plikely_R_R
    ML(1) = 2;
    ML_error = bayes_sum_R * Ppost_R_R;
else
    ML(1) = 1;
    ML_error = bayes_sum_R * Ppost_B_R;
end

if Plikely_B_B > Plikely_B_R
    ML(2) = 2;
    ML_error = ML_error + bayes_sum_B * Ppost_R_B;
else
    ML(2) = 1;
    ML_error = ML_error + bayes_sum_B * Ppost_B_B;
end

% c


% d
% Case 1:
% MAP = [2 2]
% ML = [1 2]
% MAP error: 0.1
% ML error: 0.1636

% Case 2:
% MAP = [1 2]
% ML = [1 2]
% MAP error: 0.3818
% ML error: 0.3818


%% e
% Simulate case 1 10^5 times
n = 10e5;
error_case1 = [0 0];

for i = 1:n
    error_it = simulate(1, 9, 1, 9);
    error_case1 = error_case1 + error_it;
end
error_case1 = error_case1/n;

% Estimated error: 
% MAP: 0.1003
% ML: 0.1633

% Simulate case 2 10^5 times
error_case2 = [0 0];

for i = 1:n
    error_it = simulate(4, 6, 1, 9);
    error_case2 = error_case2 + error_it;
end
error_case2 = error_case2/n;

% Estimated error:
% MAP: 0.4549
% ML: 0.3825

%% Wrap all work done in function

function outputs = RunExperiment(r1, b1, r2, b2) 
    % prior distribution - probability of getting red/blue from urn 1
    pi_R = r1/(r1+b1);
    pi_B = b1/(r1+b1);

    % likelihood fcn - given result from urn 1, prob of pulling from urn 2
    Plikely_R_B = r2/(r2+b2+1);
    Plikely_B_B = (b2+1)/(r2+b2+1);
    Plikely_R_R = (r2+1)/(r2+b2+1);
    Plikely_B_R = b2/(r2+b2+1);
    
    bayes_sum_R = Plikely_R_B*pi_B + Plikely_R_R*pi_R;
    bayes_sum_B = Plikely_B_B*pi_B + Plikely_B_R*pi_R;
    
    % posteriori distribution - given urn2, return urn1
    Ppost_R_B = Plikely_B_R*pi_R/bayes_sum_B;
    Ppost_R_R = Plikely_R_R*pi_R/bayes_sum_R;
    Ppost_B_R = Plikely_R_B*pi_B/bayes_sum_R;
    Ppost_B_B = Plikely_B_B*pi_B/bayes_sum_B;
    
    % b
    % make decision vectors
    MAP = [0 0];
    ML = [0 0];
    
    % MAP maximizes a-posteriori distribution
    if Ppost_B_R > Ppost_R_R
        MAP(1) = 2;
    else
        MAP(1) = 1;
    end
    
    if Ppost_B_B > Ppost_R_B
        MAP(2) = 2;
   else
        MAP(2) = 1;
   end
    
    % ML maximizes likelihood fcn
    if Plikely_R_B > Plikely_R_R
        ML(1) = 2;
    else
        ML(1) = 1;
   end
    
    if Plikely_B_B > Plikely_B_R
        ML(2) = 2;
    else
        ML(2) = 1;
     end

    outputs = [MAP; ML];
end

%%
function error = simulate(r1, b1, r2, b2)
    error = [0, 0];
    % choose a ball from the urn and move it
    ball = rand;
    threshold = r1/(r1+b1);
    
    if ball < threshold
        ball1 = 1; % red
        r1 = r1 - 1;
        r2 = r2 + 1;
    else
        ball1 = 2; % blue
        b1 = b1 - 1;
        b2 = b2 + 1;
    end
    
    % choose a ball from urn 2
    ball = rand;
    threshold = r2/(b2+r2);
    if ball < threshold
        ball2 = 1; % red
    else
        ball2 = 2; % blue
    end
    
    % make a decision based on MAP and ML
    decisions = RunExperiment(r1,b1,r2,b2);
    
    % check result
    if decisions(1, ball2) ~= ball1
        error(1) = error(1) + 1;
    end
    
    if decisions(2, ball2) ~= ball1
        error(2) = error(2) + 1;
    end
end
