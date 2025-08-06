%% Problem 2b: Empirical Distribution of Refractory-Constrained Spike Patterns
% Discretize 100 ms into n=10 bins, k=3 spikes, no adjacent spikes.

clear; clc;
chi2cdf = @(x,df) gammainc(x./2, df./2, 'lower');
% 1) Enumerate all valid patterns
n = 10; 
k = 3;
comb = nchoosek(1:n, k);                 % all choices of k bins
valid = comb(all(diff(comb,1,2) > 1, 2), :);  % enforce refractory (no adjacent)

numPatterns = size(valid, 1);            % should be 56

% 2) Convert to binary matrix (patterns × time bins)
patterns = zeros(numPatterns, n);
for i = 1:numPatterns
    patterns(i, valid(i,:)) = 1;
end

% 3) Sample N = 1000 patterns uniformly
N = 1000;
rep = ceil(N/numPatterns);               % how many repeats to cover N
pool = repmat(1:numPatterns, 1, rep);    % pool of indices
pool = pool(randperm(length(pool)));     % shuffle using randperm
sampleIdx = pool(1:N);                   % take first N samples

% 4) Count occurrences and normalize to get probabilities
counts = histcounts(sampleIdx, 1:(numPatterns+1));
prob = counts / N;

% 5) Plot the empirical probability distribution
figure('Name','Empirical Pattern Distribution','NumberTitle','off');
bar(prob, 'FaceColor',[0.2 0.6 0.8]);
xlabel('Pattern ID');
ylabel('Empirical Probability');
title('Uniform Sampling of Refractory-Constrained Spike Patterns (N=1000)');
ylim([0, max(prob)*1.2]);

%verify mean ≈1/numPatterns
fprintf('Mean empirical probability = %.3f (expected ≈%.3f)\n', ...
        mean(prob), 1/numPatterns);


% Extension for confirmation
% A) 95% binomial‐proportion CIs
se = sqrt( prob .* (1 - prob) / N );       % standard error for each bar
ci95 = 1.96 * se;                         % 95% z‐interval

% Plot with error bars
figure('Name','Pattern Probabilities ±95% CI','NumberTitle','off');
bar(prob, 'FaceColor',[0.2 0.6 0.8]); hold on;
errorbar(1:numPatterns, prob, ci95, 'k.', 'LineWidth',1);
hold off;
xlabel('Pattern ID');
ylabel('Empirical Probability');
title('Uniform Sampling w/ 95% CI (N=1000)');
ylim([0, max(prob+ci95)*1.2]);

% B) Chi-square goodness-of-fit test
expected_count = N/numPatterns;
chi2_stat = sum((counts - expected_count).^2 / expected_count);
df = numPatterns - 1;
p_val = 1 - chi2cdf(chi2_stat, df);

fprintf('Chi2 GOF: χ²=%.2f, df=%d, p=%.3f\n', chi2_stat, df, p_val);