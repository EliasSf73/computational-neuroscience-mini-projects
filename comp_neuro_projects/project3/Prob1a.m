% Problem 1a: simulate and plot histogram of firing rates

% paramteres
mu_plus=30;
sigma_plus=3;
mu_minus=10;
sigma_minus=3;
n_trials=100;
n_ups=n_trials/2;
n_downs=n_trials/2;

% simulate the firing rates
rng(0);
r_plus= mu_plus + sigma_plus*randn(n_ups,1); % row vector of responses
r_minus= mu_minus + sigma_minus*randn(n_downs,1); 
% combine into single vector
r=[r_plus; r_minus];
labels= [ones(n_ups,1); zeros(n_downs,1)];
%   labels==1 → true “+”, labels==0 → true “–”


% plotting the overlaid histogram
figure;
hold on;
histogram(r_plus, 'BinWidth',2, 'FaceAlpha',0.6, 'EdgeColor','none', 'FaceColor',[0 .5 0], 'DisplayName','Upward (+)');
histogram(r_minus,'BinWidth',2, 'FaceAlpha',0.6, 'EdgeColor','none', 'FaceColor',[.8 0 0], 'DisplayName','Downward (–)');
xlabel('Firing rate (Hz)');
ylabel('Count');
title('Simulated firing-rates: upwards vs downwards (100 trials)');
legend('Location','Best');
grid on;
hold off;