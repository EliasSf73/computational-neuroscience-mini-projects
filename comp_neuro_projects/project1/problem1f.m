%% Parameters
% rng(9);  % o fixed seed for reproducibility

dt = 0.1;             % time step (ms)
T = 10000;             % total simulation time (ms) = 10 second
t_array = 0:dt:T;
nSteps = length(t_array);

Erest = -65;          % resting potential (mV)
Eth = -55;            % threshold (mV)
tau = 20;             % membrane time constant (ms)
R = 10;               % membrane resistance (Ohm)
I_inj = 0;            % no constant injected current
sigma_noise = 8.5;    

nTrials = 5;          %  5 repeated trials

%% arrays for  storing results
spike_times_all = cell(nTrials, 1);  
firing_rates = zeros(nTrials, 1);    

for trial = 1:nTrials
    
    rng(trial);
    
    % Initialize membrane potential
    v = Erest * ones(1, nSteps);
    spike_times = [];
    
    % Euler integration for each time step
    for i = 1:nSteps-1
        I_noise = sigma_noise * randn;  % Gaussian noise, mean 0, std = sigma_noise
        v(i+1) = v(i) + (dt/tau)*(Erest - v(i) + R*I_noise);
        
        % threshold crossing check
        if v(i+1) >= Eth
            spike_times = [spike_times, t_array(i)];
            v(i+1) = Erest;  % reset
        end
    end
    
    %storing spike times
    spike_times_all{trial} = spike_times;
    
    % calculate firing rate for the current trial (spikes / second)
    num_spikes = length(spike_times);
    firing_rates(trial) = num_spikes / (T/1000);  % T is in ms -> T/1000 is 1 second
end

%% raster plot for all 5 trials
figure; hold on;
for trial = 1:nTrials
    these_spikes = spike_times_all{trial};
    % plotting vertical lines for each spike
    for s = 1:length(these_spikes)
        line([these_spikes(s) these_spikes(s)], [trial-0.4 trial+0.4], 'Color','k');
    end
end
xlabel('Time (ms)');
ylabel('Trial #');
title('Raster Plot: 5 Trials (T = 1 sec each)');
ylim([0 nTrials+1]);
hold off;

%% Calculate variance for firing rates 
mean_firing_rate = mean(firing_rates);
var_firing_rate = var(firing_rates);

fprintf('Firing rates for each trial: %s\n', num2str(firing_rates'));
fprintf('Mean firing rate: %.2f Hz\n', mean_firing_rate);
fprintf('Variance of firing rates: %.2f (Hz^2)\n', var_firing_rate);
