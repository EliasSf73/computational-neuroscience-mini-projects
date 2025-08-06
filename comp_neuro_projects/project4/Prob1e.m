% compare_spike_code_entropies.m
%----------------------------------------------------------------------
% Compare the Shannon entropy of Poisson vs. Exponential (geometric)
% spike-count codes as a function of the mean count λ = r·T.


clear; clc; close all;

%% 1) Define range of mean counts
% We sweep λ from 0.1 up to 100 on a logarithmic grid
lambda = logspace(-1, 2, 200);  

%% 2) Compute analytic entropy for Poisson code
% Stirling-based approximation:
%   S_pois(λ) ≈ ½ [log2(λ) + log2(2π) + log2(e)]
entropy_poisson = 0.5 * ( log2(lambda) + log2(2*pi) + log2(exp(1)) );

%% 3) Compute exact entropy for Exponential (geometric) code
%   S_exp(λ) = log2(1+λ) + λ·log2(1 + 1/λ)
entropy_exponential = log2(1 + lambda) + lambda .* log2(1 + 1 ./ lambda);

%% 4) Plot both curves on log–log axes
figure;
loglog(lambda, entropy_poisson,    'b-',  'LineWidth', 2); hold on;
loglog(lambda, entropy_exponential, 'r--', 'LineWidth', 2);
hold off;
grid on;

xlabel('\lambda = r\,T',          'FontSize', 12);
ylabel('Entropy S \rm(bits)',    'FontSize', 12);
title ('Entropy of Poisson vs. Exponential Spike-Count Codes', ...
       'FontSize', 14);

legend('Poisson code', 'Exponential code', ...
       'Location', 'northwest', 'FontSize', 12);

%% 5) Print a few example values to the console
example_lambdas = [0.5, 1, 10, 100];
fprintf('\nExample entropy values:\n');
fprintf('  λ     S_{pois} (bits)   S_{exp} (bits)\n');
fprintf('  --------------------------------------\n');
for L = example_lambdas
    S_p = 0.5 * ( log2(L) + log2(2*pi) + log2(exp(1)) );
    S_e = log2(1 + L) + L * log2(1 + 1/L);
    fprintf('  %4.1f      %7.3f            %7.3f\n', L, S_p, S_e);
end
fprintf('\n');
