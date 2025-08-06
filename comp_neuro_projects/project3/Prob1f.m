% Problem 1(f): Compare easy, moderate, and hard conditions

% parameters for each condition ---
conds = { ...
  struct('name','easy',     'mu_p',30,'sigma_p',3, 'mu_m',10,'sigma_m',3), ...
  struct('name','moderate','mu_p',26,'sigma_p',4, 'mu_m',15,'sigma_m',4), ...
  struct('name','hard',    'mu_p',20,'sigma_p',5, 'mu_m',18,'sigma_m',5)  ...
};

n_trials = 100;     % total per condition
n_up     = n_trials/2;
n_down   = n_up;
Z        = 0:1:50;  % fine‐grained thresholds

%storage for plotting & optimal values
colors = [ [0 .7 0]; [0 .4 .9]; [1 .2 0] ];  % green, blue, red
figure; hold on;
optimal = struct();

% Loop over conditions 
for c = 1:numel(conds)
  p = conds{c};
  
  % a) Simulate firing rates
  rng(0);  % keep same seed for fair comparison
  r_up   = p.mu_p    + p.sigma_p*randn(n_up,1);
  r_down = p.mu_m    + p.sigma_m*randn(n_down,1);
  r      = [r_up; r_down];
  labels = [ones(n_up,1); zeros(n_down,1)];  
  
  % b) Sweep thresholds
  beta_vals  = zeros(size(Z));
  alpha_vals = zeros(size(Z));
  p_vals     = zeros(size(Z));
  for k = 1:length(Z)
    dec = (r > Z(k));
    beta_vals(k)  = sum(dec(labels==1)) / n_up;
    alpha_vals(k) = sum(dec(labels==0)) / n_down;
    p_vals(k)     = (beta_vals(k) + (1-alpha_vals(k)))/2;
  end

  % c) Plot ROC
  plot(alpha_vals, beta_vals, '-o', 'Color', colors(c,:), ...
       'DisplayName', sprintf('%s (%.0f/%.0f)', p.name, p.mu_p, p.mu_m));
  
  % d) Find optimal threshold & accuracy
  [p_best, idx]         = max(p_vals);
  z_best                = Z(idx);
  optimal(c).name       = p.name;
  optimal(c).z_opt      = z_best;
  optimal(c).accuracy   = p_best;
end

%ROC plot ---
plot([0 1],[0 1],'k--','DisplayName','chance');
xlabel('False‐alarm rate \alpha');
ylabel('Hit rate \beta');
title('ROC Curves: Easy, Moderate, and Hard Conditions');
legend('Location','SouthEast');
axis([0 1 0 1]);
grid on;
