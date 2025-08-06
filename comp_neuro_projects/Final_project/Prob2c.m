%% Problem 2c: Entropy Calculation S1 [bits]
% Discretize 100 ms into 10 bins, exactly 3 spikes, 10 ms refractory.
% Uniformly sample N=1000 patterns, compute empirical probabilities, then entropy.

clear; close all; clc;

%% 1) Enumerate all valid refractory‐constrained patterns
n = 10;    % number of bins
k = 3;     % spikes
allComb = nchoosek(1:n, k);                          % all ways to pick k bins
validIdx = all( diff(allComb,1,2) > 1, 2 );          % reject adjacent picks
validPatterns = allComb(validIdx, :);                % size P×k
P = size(validPatterns, 1);                          % should be 56

%% 2) Convert to binary matrix (P×n)
patterns = zeros(P, n);
for i = 1:P
    patterns(i, validPatterns(i,:)) = 1;
end

%% 3) Draw N=1000 samples uniformly
N = 1000;
reps = ceil(N/P);
pool = repmat(1:P, 1, reps);
pool = pool(randperm(numel(pool)));    % shuffle
sampleIdx = pool(1:N);

%% 4) Compute empirical probabilities
counts = histcounts(sampleIdx, 1:(P+1));
prob    = counts / N;                  % 1×P

%% 5) Compute Shannon entropy S1 [bits]
% ignore zero-prob entries (should be none here)
nonzero = prob > 0;
S1 = -sum( prob(nonzero) .* log2( prob(nonzero) ) );

fprintf('Empirical entropy S1 = %.4f bits\n', S1);
fprintf('Theoretical max log2(P) = %.4f bits\n', log2(P));

%% 6)  Visualize the distribution
figure('Name','Empirical Pattern Distribution','NumberTitle','off');
bar(prob, 'FaceColor',[0.2 0.6 0.8]);
xlabel('Pattern ID');
ylabel('Empirical Probability');
title(sprintf('Uniform Sampling (N=%d), Entropy S1=%.4f bits', N, S1));
ylim([0, max(prob)*1.2]);
