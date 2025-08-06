

% -------------------------------------------------------------------------
%  HW#4  Problem 1(a)
%  Simulate Poisson spike counts and plot empirical P(n)
% -------------------------------------------------------------------------


% Parameters
r=15; % firing rate for each trial (Hz)
T=1; % window length (S)
N=100; % number of trials
lambda= r*T; % rate overall trial time

% Let X~ the number of spikes for given trial: P(x=n) is poisson
% Let's draw poisson counts for N trials
rng(41);
counts= poissrnd(lambda, [N,1]); % 100 poisson draws

% converting counts to an empirical pmf

n_max   = max(counts);          % highest count observed
edges   = -0.5 : 1 : (n_max+0.5);           % unit-width bin edges
[count_hist, ~] = histcounts(counts, edges);% bin the data and count spikes within each bin
p_emp   = count_hist / N;       % normalise counts in each bin to probabilities
n_vals  = 0 : n_max;            % n corresponding to p_emp(k)

% Histogram Plotting
figure;
bar(n_vals, p_emp, 'FaceColor', [0.2 0.5 0.8]);
xlabel('Spike count  n');
ylabel('Empirical  P(n)');
title(sprintf('Poisson spike counts: r = %g Hz, T = %gs, N = %d', r, T, N));
grid on;

% basic stats
fprintf('Mean of simulated counts  : %.2f spikes\n', mean(counts));
fprintf('Theoretical Poisson mean  : %.2f spikes\n', lambda);