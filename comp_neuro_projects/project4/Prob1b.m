
% -------------------------------------------------------------------------
%  HW#4  Problem 1(b)
%  Entropy Approximation: Numerical Vs analytical entropy
% -------------------------------------------------------------------------

% GETTING N_VALS AND P_EMP FROM PROBLEM 1(a)
%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------

% Numerical entropy 
  idx    = p_emp > 0; % (ignore zero-probability bins)
  S_num  = -sum( p_emp(idx) .* log2(p_emp(idx)) );

 % Estimate mean-count(mc) from the pmf
 mc= sum(n_vals.*p_emp);
 S_an= 0.5* (log2(mc)  + log2(2*pi)   + log2(exp(1))  );

 % absolute error between the two
 delta= abs(S_num-S_an);

 % Display
  fprintf('Numerical entropy: %.3f bits\n', S_num);
  fprintf('Analytic entropy : %.3f bits\n', S_an);
  fprintf('Absolute error    : %.3f bits\n', delta);
 