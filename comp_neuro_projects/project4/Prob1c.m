
% Prob1c_sweep.m
% Study how the absolute entropy error ΔS scales with
%  1) number of trials N (at fixed λ)
%  2) mean count λ (at fixed N)
% and fit power‐law slopes in the appropriate regimes.

clear; close all; clc;


% --- Part A: ΔS vs N, fixed lambda = 15 ---
r = 15;      % Hz
T = 1;       % s
lambda = r*T;
N_vals = round(logspace(2,5,10));  % from 1e2 to 1e5
deltaN = zeros(size(N_vals));

for i = 1:numel(N_vals)
    rng(0);
    N = N_vals(i);
    % simulate and get pmf
    counts = poissrnd(lambda, [N,1]);
    edges  = -0.5:1:(max(counts)+0.5);
    phist  = histcounts(counts, edges)/N;
    n_vals = 0:(numel(phist)-1);
    % numerical entropy
    idx    = phist>0;
    S_num  = -sum(phist(idx).*log2(phist(idx)));
    % analytic entropy
    S_an   = 0.5*(log2(lambda)+log2(2*pi)+log2(exp(1)));
    deltaN(i) = abs(S_num - S_an);
end

figure;
loglog(N_vals, deltaN, '-o','LineWidth',1.5);
xlabel('Number of trials  N');
ylabel('\Delta S (bits)');
title('\DeltaS vs N  ( \lambda = 15 )');
grid on;

%  fit a line to check slope ~ -1/2
p = polyfit(log10(N_vals), log10(deltaN), 1);
fprintf('Slope vs N (log-log fit): %.2f (should be ≈ -0.5)\n', p(1));


% --- Part B: ΔS vs lambda, fixed N = 10000 ---
N = 1e4;
T = 1;  % second
lambda_vals = [5 10 20 50 100];  % try a range
deltaL = zeros(size(lambda_vals));

for i = 1:numel(lambda_vals)
    lambda = lambda_vals(i);
    r = lambda / T;
    counts = poissrnd(lambda, [N,1]);
    edges  = -0.5:1:(max(counts)+0.5);
    phist  = histcounts(counts, edges)/N;
    n_vals = 0:(numel(phist)-1);

    idx    = phist>0;
    S_num  = -sum(phist(idx).*log2(phist(idx)));
    S_an   = 0.5*(log2(lambda)+log2(2*pi)+log2(exp(1)));
    deltaL(i) = abs(S_num - S_an);
end

figure;
loglog(lambda_vals, deltaL, '-s','LineWidth',1.5);
xlabel('\lambda = rT');
ylabel('\Delta S (bits)');
title('\DeltaS vs \lambda  (N = 10^4)');
grid on;

% Optional: fit slope vs lambda
q = polyfit(log10(lambda_vals), log10(deltaL), 1);
fprintf('Slope vs lambda (log-log fit): %.2f (should be ≈ -1)\n', q(1));
