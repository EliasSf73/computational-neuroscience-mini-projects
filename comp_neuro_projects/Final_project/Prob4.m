%% Revised HW3–1b (toolbox-free) : Analytic curves + stochastic demo
%% Go / No-go decoding : original "easy" parameters
clear; close all; clc;

% -------- parameters (unchanged from HW3) -------------------------------
muP = 30;          % mean rate for "+"  (Hz)
muM = 10;          % mean rate for "–"  (Hz)
sig =  3;          % common std-dev     (Hz)

% -------- helper: Gaussian CDF without toolbox --------------------------
Phi = @(x) 0.5*(1 + erf(x./sqrt(2)));

% -------- analytic curves -----------------------------------------------
z = 0:0.2:40;                             % threshold axis
beta  = 1 -   Phi((z-muP)/sig);           % hit rate
alpha = 1 -   Phi((z-muM)/sig);           % false-alarm rate
p     = 0.5*(beta + 1 - alpha);           % overall accuracy

[~,idx] = max(p);           % analytic optimum -> 20 Hz exactly
z_opt  = z(idx);            % should equal 20

% -------- stochastic 100-trial demo (random ± mix) ----------------------
N = 100;
isPlus = rand(1,N) < 0.5;                 % Bernoulli(0.5) each trial
rates  = isPlus .* (muP + sig*randn(1,N)) + ...
         ~isPlus.* (muM + sig*randn(1,N));

z_fixed = 20;                             % fixed HW3 threshold
hits  = sum(rates(isPlus)  > z_fixed);
fas   = sum(rates(~isPlus) > z_fixed);
beta_emp  = hits  / max(1,sum(isPlus));    % avoid /0 if class is empty
alpha_emp = fas   / max(1,sum(~isPlus));
p_emp     = 0.5*(beta_emp + 1 - alpha_emp);

fprintf('Empirical (random ±, z=20):  beta=%.2f  alpha=%.2f  p=%.2f\n',...
        beta_emp, alpha_emp, p_emp);

% -------- plot ----------------------------------------------------------
figure; hold on;
plot(z, alpha, 'r-', 'DisplayName','false alarm  \alpha(z)');
plot(z, beta,  'b-', 'DisplayName','hit rate      \beta(z)');
plot(z, p,     'k-', 'LineWidth',1.4, 'DisplayName','accuracy     p(z)');
plot(z_opt, p(idx), 'ko', 'MarkerFaceColor','g', 'DisplayName','optimal z^*');
xlabel('threshold  z  (Hz)'); ylabel('probability');
legend('Location','best'); grid on;
title('Go / No-go curves  (\mu_+=30, \mu_-=10, \sigma=3)');

