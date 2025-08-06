%% Problem 3c: Cost Estimation for Random (a,b) Choices
% Draw several random decision boundaries r2 = a*r1 + b, classify the
% four AND inputs, and measure the total misclassification cost E.

clear; close all; clc;

% 1) AND inputs and desired outputs
inputs  = [0 0;  0 1;  1 0;  1 1];  % [r1, r2]
desired = all(inputs, 2);           % AND truth: [0;0;0;1]

% 2) Randomly sample K (a,b) pairs
K      = 5;                       
a_vals = -3 + 6*rand(K,1);         % slopes in [-3,3]
b_vals = -1 + 4*rand(K,1);         % intercepts in [-1,3]

% 3) Compute error E(a,b) = # of misclassified points
fprintf(' trial   a      b     cost E\n');
fprintf('-----------------------------\n');
for j = 1:K
    a = a_vals(j);
    b = b_vals(j);
    
    % classify by line: r_out=1 if r2>a*r1+b
    r_out = inputs(:,2) > (a*inputs(:,1) + b);
    
    % total absolute difference from desired AND outputs
    E = sum(abs(r_out - desired));
    
    fprintf('  %1d     %+5.2f  %+5.2f    %d\n', j, a, b, E);
end


%% Extension: Classification Margins for the AND Boundary
% For the true AND line (a=-1,b=1.5), compute how far each input point
% lies from the decision boundary in normalized distance units.

% 1) Reconstruct perceptron parameters
w1    =  1;      % from a=-1 ⇒ w1=-a*w2 with w2=1
w2    =  1;
Theta =  1.5;    

% 2) Convert AND labels to ±1
X = inputs;                
y = 2*desired - 1;         % 0→-1, 1→+1

% 3) Raw activations u = w·x - Θ
u = X(:,1)*w1 + X(:,2)*w2 - Theta;

% 4) Signed distance to boundary: d = u / ||w||
norm_w = sqrt(w1^2 + w2^2);
d      = u / norm_w;

% 5) Report margins
fprintf('\n r1 r2  label  activation   margin d\n');
fprintf('--------------------------------------\n');
for i = 1:4
    fprintf(' %d  %d   %2d     %+6.2f     %+6.2f\n', ...
            X(i,1), X(i,2), y(i), u(i), d(i));
end
min_margin = min(abs(d));
fprintf('\nMinimum classification margin = %.3f units\n', min_margin);

%% Error vs. Random (a,b), and Classification Margins Bar Plot
%% 2) Sample K random (a,b) and compute cost E
K      = 50;                         % increase samples for a smoother scatter
a_vals = -3 + 6*rand(K,1);          % slope in [-3,3]
b_vals = -1 + 4*rand(K,1);          % intercept in [-1,3]
E_vals = zeros(K,1);

for j = 1:K
    a = a_vals(j);
    b = b_vals(j);
    % classify via line: r_out=1 if r2>a*r1+b
    r_out    = inputs(:,2) > (a*inputs(:,1) + b);
    % total misclassifications
    E_vals(j)= sum(abs(r_out - desired));
end

%% 3) Scatter plot of (a,b) colored by error E
figure('Name','3c: Error vs. (a,b)','NumberTitle','off');
scatter(a_vals, b_vals, 100, E_vals, 'filled');
colorbar;
colormap(parula);
caxis([0,4]);
xlabel('slope a');
ylabel('intercept b');
title('Misclassification Error E(a,b) for Random Lines');
grid on;

%% Extension Plot: Classification Margins for the true AND boundary
% Reconstruct perceptron params from a=-1,b=1.5
w1 = 1;  w2 = 1;  Theta = 1.5;

% Compute signed distances d_i for each AND input
u = inputs(:,1)*w1 + inputs(:,2)*w2 - Theta;  % raw activation
norm_w = sqrt(w1^2 + w2^2);
d = u ./ norm_w;                               % signed margin

%% 4) Bar plot of margins
figure('Name','3c Extension: Classification Margins','NumberTitle','off');
bar(d, 0.6, 'FaceColor',[0.2 0.6 0.8]);
hold on;
yline(0, 'k--', 'LineWidth',1);   % decision threshold
hold off;
xticks(1:4);
xticklabels({'(0,0)','(0,1)','(1,0)','(1,1)'});
xlabel('Input pair (r_1,r_2)');
ylabel('Signed distance to boundary');
title('Classification Margins for AND Boundary');
grid on;