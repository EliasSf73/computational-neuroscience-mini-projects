%% Problem 2f: Empirical Distribution under Poisson Sampling
% Discretize a 100 ms window into n=10 bins, with Poisson rate 30 Hz ⇒
% p = λ·Δt = 30 Hz×0.010 s = 0.3 chance of ≥1 spike in each bin.
% We draw N=1000 independent patterns and plot their empirical probabilities.

clear; close all; clc;

%% 1) Parameters
n = 10;              % number of bins (100 ms/10 ms)
N = 1000;            % number of patterns to sample
p = 0.3;             % spike probability per bin under Poisson

%% 2) Generate spike patterns (N×n binary matrix)
spikes = rand(N, n) < p;  

%% 3) Encode each pattern as an integer ID from 1 to 2^n
%    ID = 1 + sum_{j=1..n} spikes(i,j)*2^(n-j)
weights = 2.^(n-1:-1:0);        
ids = spikes * weights.' + 1;  

%% 4) Compute empirical probability of each of the 2^n patterns
P = 2^n;                           
counts = histcounts(ids, 1:(P+1));
prob_emp = counts / N;            

%% 5) Plot the empirical distribution
figure('Name','Poisson Pattern Distribution','NumberTitle','off');
bar(prob_emp, 'FaceColor',[0.3 0.7 0.4]);
xlabel('Pattern ID');
ylabel('Empirical Probability');
title(sprintf('Poisson Sampling (p=%.2f), N=%d', p, N));
xlim([1 P]);
ylim([0, max(prob_emp)*1.2]);

%% extension: overlay theoretical pattern probabilities
% theoretical P(pattern) = p^k * (1-p)^(n-k), where k = # of spikes in that pattern
theo_prob = zeros(1,P);
for id = 1:P
    bits = bitget(id-1, n:-1:1);
    k = sum(bits);
    theo_prob(id) = p^k * (1-p)^(n-k);
end
hold on;
plot(1:P, theo_prob, '.r', 'MarkerSize',8);
hold off;
legend('Empirical','Theoretical','Location','northeast');

%%%Extension: Spike‐Count Distribution vs. Binomial Theory (fixed)

% 'spikes' is our N×n matrix from above
k_counts = sum(spikes, 2);                          % total spikes per trial
counts_k  = histcounts(k_counts, -0.5:1:10.5);       % bins centered 0,1,…,10
emp_k     = counts_k / N;                            % empirical P(k)

% Theoretical Binomial(10,0.3) via factorials
n_bins = 10; p = 0.3;
k_vals = 0:n_bins;
% combination n choose k for a vector:
comb = factorial(n_bins) ./ (factorial(k_vals) .* factorial(n_bins - k_vals));
theo_k = comb .* (p.^k_vals) .* ((1-p).^(n_bins - k_vals));

% Plot
figure('Name','Spike‐Count Distribution','NumberTitle','off');
bar(k_vals, emp_k, 0.6, 'FaceColor',[0.2 0.6 0.8]); hold on;
plot(k_vals, theo_k, 'r-o', 'LineWidth',1.5, 'MarkerSize',6);
hold off;
xlabel('Number of spikes in 10 bins (k)');
ylabel('Probability');
title('Empirical vs. Theoretical Binomial(10,0.3) Distribution');
legend('Empirical','Binomial Theory','Location','Best');
grid on;
