

clear; clc;

% parameters
N      = 10;          % trials
T      = 0.1;        % s
dt     = 1e-3;        % s  (1 ms bin)
nbins  = round(T/dt); % 100 bins
rate_list = 5:5:50;   % 5,10,…,50 Hz
num_rates = numel(rate_list);

mean_rate = zeros(1,num_rates);
var_rate  = zeros(1,num_rates);

rng(1);              % repeatability

for i = 1:num_rates
    r = rate_list(i);           % <-- correct indexing
    p = r*dt;                   % Bernoulli prob. per bin
    spk = rand(nbins,N) < p;    % logical spikes (rand not randn)

    trial_rate     = sum(spk,1)/T;   % spikes / 0.1 s  → Hz
    mean_rate(i)   = mean(trial_rate);
    var_rate(i)    = var(trial_rate,1); % population variance (flag = 1)
end

% log‑log scatter
figure;  loglog(mean_rate, var_rate, 'ko', 'MarkerFaceColor','k'); hold on;
xlabel('Mean firing rate \mu (Hz)');
ylabel('Variance \sigma^2 (Hz^2)');
grid on; box on;

% linear fit in log–log space
p        = polyfit(log(mean_rate), log(var_rate), 1);
slope    = p(1);

intercept= p(2);
Fano     = exp(intercept);                      % intercept = log F
loglog(mean_rate, exp(intercept)*mean_rate.^slope, 'r-', 'LineWidth',1.5);
legend('data','log–log fit','Location','NorthWest');
title(sprintf('Fano factor estimate  F ≈ %.3f   (slope = %.2f)', Fano, slope));

fprintf('Estimated Fano factor: %.3f\n', Fano);








