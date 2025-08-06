
% Problem 1(d): Extension of Problem 1(a) with ROC

% Parameters
mu_plus    = 30;    sigma_plus  = 3;
mu_minus   = 10;    sigma_minus = 3;
n_trials   = 10;
n_up       = n_trials/2;  % 5
n_down     = n_up;         % 5

% Simulate 10 trials
rng(0);
r_plus     = mu_plus  + sigma_plus*randn(n_up,1);
r_minus    = mu_minus + sigma_minus*randn(n_down,1);

% Combine & label
r          = [r_plus; r_minus];
labels     = [ones(n_up,1); zeros(n_down,1)];  % 1=up, 0=down

% Thresholds to test
Z          = 0:3:45;

% Preallocate arrays
beta_vals  = zeros(size(Z));
alpha_vals = zeros(size(Z));
p_vals     = zeros(size(Z));

% Steps 2–3: sweep thresholds
for k = 1:length(Z)
    z_k           = Z(k);
    decisions_k   = (r > z_k);

    % Hit rate β(z)
    hits          = sum(decisions_k(labels==1));
    beta_vals(k)  = hits / n_up;

    % False‐alarm rate α(z)
    fas           = sum(decisions_k(labels==0));
    alpha_vals(k) = fas / n_down;

    % Percent‐correct p(z)
    p_vals(k)     = (beta_vals(k) + (1 - alpha_vals(k))) / 2;
end

% Display numeric results
fprintf('Thresholds:            %s\n', mat2str(Z));
fprintf('Hit rates (β):         %s\n', mat2str(beta_vals,3));
fprintf('False‐alarm rates (α): %s\n', mat2str(alpha_vals,3));
fprintf('Percent correct (p):   %s\n', mat2str(p_vals,3));

% Step 4: Plot ROC curve (α vs β)
figure;
plot(alpha_vals, beta_vals, 'o-', 'LineWidth',1.5, 'MarkerSize',8);
hold on;
plot([0 1],[0 1],'k--');           % chance line
xlabel('False‐alarm rate \alpha');
ylabel('Hit rate \beta');
title('ROC Curve for 10‐Trial Simulation');
grid on;
axis([0 1 0 1]);
legend('Data','Chance','Location','SouthEast');

% Step 5: Plot percent‐correct vs threshold
figure;
plot(Z, p_vals, 'b-^', 'LineWidth',1.5, 'MarkerSize',8);
xlabel('Threshold z (Hz)');
ylabel('Percent correct p');
title('Percent Correct vs. Decision Threshold');
grid on;


% % Problem 1d: Extension of problem 1a
% 
% % paramteres
% mu_plus=30;
% sigma_plus=3;
% mu_minus=10;
% sigma_minus=3;
% n_trials=10;
% n_ups=n_trials/2;
% n_downs=n_trials/2;
% 
% % simulate the firing rates
% rng(0);
% r_plus= mu_plus + sigma_plus*randn(n_ups,1); % row vector of responses
% r_minus= mu_minus + sigma_minus*randn(n_downs,1); 
% % combine into single vector
% r=[r_plus; r_minus];
% labels= [ones(n_ups,1); zeros(n_downs,1)];
% %   labels==1 → true “+”, labels==0 → true “–”
% 
% % thresholds
% z=0:3:45;
% beta_vals= zeros(size(z));
% alpha_vals=zeros(size(z));
% p_vals=zeros(size(z));
% 
% for k=1:length(z)
%     z_k=z(k);
%     % decision at this threshold
%     decision_k=(r>z_k);
%     % hits rate
%     hits_k=sum(decision_k(labels==1));
%     beta_vals(k)=hits_k/n_up;
%     % false alarm rate
%     falseAlarms_k= sum(decision_k(labels==0));
%     alpha_vals(k)=falseAlarms_k/n_down;
%     % percent correct p(z)
%     p_vals(k)= (beta_vals(k) + 1-alpha_vals(k))/2;
% end
% 
% fprintf('Thresholds:           %s\n', mat2str(z));
% fprintf('Hit rates (β):        %s\n', mat2str(beta_vals,3));
% fprintf('False-alarm rates (α):%s\n', mat2str(alpha_vals,3));
% fprintf('Percent correct (p):  %s\n', mat2str(p_vals,3));
% 
% 
% figure;
% plot(alpha_vals, beta_vals, 'o-', 'LineWidth',1.5, 'MarkerSize',8);
% hold on;
% plot([0 1],[0 1],'k--');           % chance line
% xlabel('False‐alarm rate \alpha');
% ylabel('Hit rate \beta');
% title('ROC Curve for 10‐Trial Simulation');
% grid on;
% axis([0 1 0 1]);
% legend('Data','Chance','Location','SouthEast');
% 
% % Step 5: Plot percent‐correct vs threshold
% figure;
% plot(z, p_vals, 'b-^', 'LineWidth',1.5, 'MarkerSize',8);
% xlabel('Threshold z (Hz)');
% ylabel('Percent correct p');
% title('Percent Correct vs. Decision Threshold');
% grid on;