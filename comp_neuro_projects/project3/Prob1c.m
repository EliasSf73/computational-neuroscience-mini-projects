% Problem 1(c): Sweep threshold and compute α, β, p

% paramteres
mu_plus=30;
sigma_plus=3;
mu_minus=10;
sigma_minus=3;
n_trials=100;
n_ups=n_trials/2;
n_downs=n_trials-n_ups;


% simulate the firing rates
rng(0);
r_plus= mu_plus + sigma_plus*randn(n_ups,1); % row vector of responses
r_minus= mu_minus + sigma_minus*randn(n_downs,1); 
% combine into single vector
r=[r_plus; r_minus];
labels= [ones(n_ups,1); zeros(n_downs,1)];
n_up=sum(labels==1);
n_down=sum(labels==0);

% thresholds to test
z=[20,15,10,5,1];
beta_vals= zeros(size(z));
alpha_vals=zeros(size(z));
p_vals=zeros(size(z));

for k=1:length(z)
    z_k=z(k);
    % decision at this threshold
    decision_k=(r>z_k);
    % hits rate
    hits_k=sum(decision_k(labels==1));
    beta_vals(k)=hits_k/n_up;
    % false alarm rate
    falseAlarms_k= sum(decision_k(labels==0));
    alpha_vals(k)=falseAlarms_k/n_down;
    % percent correct p(z)
    p_vals(k)= (beta_vals(k) + 1-alpha_vals(k))/2;
end

fprintf('Thresholds:           %s\n', mat2str(z));
fprintf('Hit rates (β):        %s\n', mat2str(beta_vals,3));
fprintf('False-alarm rates (α):%s\n', mat2str(alpha_vals,3));
fprintf('Percent correct (p):  %s\n', mat2str(p_vals,3));

figure;
plot(z, beta_vals, 'g-o', 'LineWidth',1.5, 'DisplayName','Hit rate \beta');
hold on;
plot(z, alpha_vals,'r-s', 'LineWidth',1.5, 'DisplayName','False‐alarm \alpha');
plot(z, p_vals,    'b-^', 'LineWidth',1.5, 'DisplayName','Percent correct p');
xlabel('Threshold z (Hz)');
ylabel('Rate');
legend('Location','Best');
title('Rates vs. Decision Threshold');
grid on;