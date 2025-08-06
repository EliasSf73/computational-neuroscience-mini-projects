%% Problem 1a: Fast-Spiking ≈50 Hz Parameter Search
% Use the Izhikevich model (Euler, dt=1 ms) at I=10 mA
% to find (a,b,c,d) that yield ~50 Hz firing in fast-spiking mode.
% We run a vectorized grid search over plausible values, then
% re-simulate & plot V(t) and the ISI histogram for the best set.

clear; close all;

%% 1) Simulation settings
dt = 1;              % time step (ms)
T  = 500;            % total duration (ms)
t  = 0:dt:T;         % time vector
N  = numel(t);       
I  = 10;             % constant input current (mA)

%% 2) Parameter grid for (a,b,c,d)
a_list = [0.05, 0.10, 0.15];
b_list = [0.15, 0.20, 0.25];
c_list = -75:5:-65;
d_list =  1:4;
[AA,BB,CC,DD] = ndgrid(a_list, b_list, c_list, d_list);
P = numel(AA);                         % total candidates
params = [AA(:), BB(:), CC(:), DD(:)]; % flatten to P×4

%% 3) Initialize state for all candidates
V = -65 * ones(P,1);         % membrane voltage
u = params(:,2) .* V;        % recovery variable = b * V
spike_count = zeros(P,1);    % spike counter

%% 4) Single time-loop: update all P models in parallel
for k = 2:N
  % a) Euler update for V and u
  dV = 0.04*V.^2 + 5*V + 140 - u + I;
  du = params(:,1) .* (params(:,2).*V - u);
  V  = V + dt*dV;
  u  = u + dt*du;
  
  % b) Detect spikes and tally
  spk = (V >= 30);
  spike_count(spk) = spike_count(spk) + 1;
  
  % c) Reset spiked models
  V(spk) = params(spk,3);        % reset to c
  u(spk) = u(spk) + params(spk,4); % add d
end

%% 5) Compute firing rates and pick best ~50 Hz
rates = spike_count / (T/1000);      % convert to Hz
[~, idx] = min(abs(rates - 50));     % find closest to 50 Hz
best = params(idx,:);
f_opt = rates(idx);

fprintf('Best fast-spiking params → a=%.2f, b=%.2f, c=%d, d=%.1f  (f≈%.1f Hz)\n', ...
        best(1), best(2), best(3), best(4), f_opt);

%% 6) Re-simulate with optimal (a,b,c,d) and plot
a = best(1);  b = best(2);  c = best(3);  d = best(4);
V = -65*ones(1,N);        % reset for single model
u = b * V;                
spikes = [];

for k = 2:N
  dV = 0.04*V(k-1)^2 + 5*V(k-1) + 140 - u(k-1) + I;
  du = a*(b*V(k-1) - u(k-1));
  V(k) = V(k-1) + dt*dV;
  u(k) = u(k-1) + dt*du;
  
  if V(k) >= 30
    spikes(end+1) = t(k);   %#ok<AGROW>
    V(k) = c;               
    u(k) = u(k) + d;        
  end
end

ISIs = diff(spikes);

figure('Name','Fast-Spiking ≈50 Hz','NumberTitle','off');
subplot(2,1,1);
plot(t, V, 'LineWidth',1);
xlabel('Time (ms)');
ylabel('V (mV)');
title(sprintf('Fast-Spiking ≈%.1f Hz  (a=%.2f,b=%.2f,c=%d,d=%.1f)', ...
      f_opt, a, b, c, d));

subplot(2,1,2);
histogram(ISIs,10);
xlabel('ISI (ms)');
ylabel('Count');
title('ISI Histogram');
