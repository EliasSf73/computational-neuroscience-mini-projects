%% Problem 3a: Perceptron “AND” Gate
% Choose weights (w1,w2) and threshold Θ so that the perceptron
% implements the Boolean AND on inputs r1,r2 ∈ {0,1}.

clear; close all; clc;

% 1) Pick parameters by hand
w1    = 1;      % weight on input r1
w2    = 1;      % weight on input r2
Theta = 1.5;    % threshold

% 2) Define all four input pairs [r1 r2]
inputs = [0 0;
          0 1;
          1 0;
          1 1];

% 3) Compute the perceptron output for each
net_input = inputs(:,1)*w1 + inputs(:,2)*w2;   % weighted sum
r_out     = net_input > Theta;                 % Boolean step

% 4) Display results alongside the AND truth table
disp(' r1  r2   net_input   perceptron_out   AND_expected')
for i = 1:size(inputs,1)
    r1 = inputs(i,1);
    r2 = inputs(i,2);
    actual = r_out(i);
    expected = all([r1 r2]);  % AND: 1 only if both r1 and r2 are 1
    fprintf('  %d    %d       %.1f           %d              %d\n', ...
            r1, r2, net_input(i), actual, expected);
end
