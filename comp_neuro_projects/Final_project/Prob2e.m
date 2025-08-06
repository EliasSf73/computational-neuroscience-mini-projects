%% Problem 2e: Standalone Pattern‐Count for Unconstrained Poisson Case
% Compute the total number of possible 10‐bin spike patterns when
% each bin can independently be 0 or 1 (no refractory or fixed‐count).

clear; clc;

n = 10;                % number of 10 ms bins in 100 ms
P = 2^n;               % total binary patterns

fprintf('Total possible patterns (2^%d) = %d\n', n, P);
