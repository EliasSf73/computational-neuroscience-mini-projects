% Problem 1b: Hit rate and false-alarm rate

% paramteres
mu_plus=30;
sigma_plus=3;
mu_minus=10;
sigma_minus=3;
n_trials=100;
n_ups=n_trials/2;
n_downs=n_trials-n_ups;
z=20;

% simulate the firing rates
rng(0);
r_plus= mu_plus + sigma_plus*randn(n_ups,1); % row vector of responses
r_minus= mu_minus + sigma_minus*randn(n_downs,1); 
% combine into single vector
r=[r_plus; r_minus];
labels= [ones(n_ups,1); zeros(n_downs,1)];
%   labels==1 → true “+”, labels==0 → true “–”

% Decision Threshold
decisions= (r>z);   % 1 = decide “+”, 0 = decide “–”

% Hit Rate
n_up = sum(labels==1); %  number of true “+” trials
hits= sum(decisions(labels==1)); % how many of those we marked 't'
beta=hits/n_up; % hit rate

% False alarm rate
n_down= sum(labels==0); % number of true '-' trials

falseAlarms=sum(decisions(labels==0)); % how many '-' ves we wrongly called '+'
alpha= falseAlarms/n_down;

% P: overall correct percent
p_correct= (beta + (1-alpha))/2;

fprintf('Hit rate (β):         %.2f\n', beta);
fprintf('False‐alarm rate (α):  %.2f\n', alpha);
fprintf('Percent correct (p):   %.2f\n', p_correct);

 





