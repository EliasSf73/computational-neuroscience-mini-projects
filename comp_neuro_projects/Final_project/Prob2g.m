%% Problem 2g: Standalone Entropy Calculation S2 [bits]
% Using the empirical distribution from Poisson sampling (Problem 2f),
% compute the Shannon entropy S2 in bits.

clear; clc;

% 1) Poisson sampling parameters & regenerate patterns
n = 10;             % number of bins
N = 1000;           % number of patterns
p = 0.3;            % spike probability per bin

% 2) Generate N Poisson‐sampled patterns
spikes = rand(N, n) < p;

% 3) Encode patterns as IDs 1..2^n
weights = 2.^(n-1:-1:0);
ids = spikes * weights.' + 1;

% 4) Empirical distribution over all 2^n patterns
P = 2^n;
counts = histcounts(ids, 1:(P+1));
prob_emp = counts / N;

% 5) Compute Shannon entropy S2 [bits]
nonzero = prob_emp > 0;
S2 = -sum(prob_emp(nonzero) .* log2(prob_emp(nonzero)));

fprintf('Entropy S2 (Poisson sampling) = %.4f bits\n', S2);
