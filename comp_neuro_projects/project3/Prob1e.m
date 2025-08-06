% Monte Carlo for required n
mu_plus   = 30; sigma_plus  = 3;
mu_minus  = 10; sigma_minus = 3;
delta     = 0.5;       % desired accuracy in Hz
alphaCI   = 0.05;      % 95% confidence
z_true    = (mu_plus+mu_minus)/2;

% candidate sample sizes per condition
candidates = 5:5:200;  % try n=5,10,...,200
nIter      = 5000;     % Monte Carlo iterations

results = zeros(size(candidates));

for i = 1:length(candidates)
  n = candidates(i);
  countHit = 0;
  for it = 1:nIter
    % draw n up‐trials and n down‐trials
    up = mu_plus  + sigma_plus*randn(n,1);
    dn = mu_minus + sigma_minus*randn(n,1);
    % estimate z_hat = midpoint of sample means
    z_hat = (mean(up)+mean(dn))/2;
    if abs(z_hat - z_true) <= delta
      countHit = countHit + 1;
    end
  end
  results(i) = countHit / nIter;  % empirical probability
end

% find smallest n with ≥95% success
idx = find(results >= 1 - alphaCI, 1, 'first');
fprintf('Required n per direction ≈ %d ( %.2f reliability )\n', ...
        candidates(idx), results(idx));

% plot reliability vs n
figure;
plot(candidates, results, 'b-o', 'LineWidth',1.5);
hold on;
yline(1-alphaCI,'r--','95% target');
xlabel('Trials per condition (n)');
ylabel('P(|\hat z - z^*|\le \delta)');
title(sprintf('Monte Carlo for \\delta=%.2f Hz, 95%% CI',delta));
grid on;
