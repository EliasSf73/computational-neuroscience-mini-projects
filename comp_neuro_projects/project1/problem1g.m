% define parameters
% rng(9);  % o fixed seed for reproducibility

dt = 0.1;             % time step (ms)
T = 3000;             % total simulation time (ms) = 1 second
t_array = 0:dt:T;
nSteps = length(t_array);

Erest = -65;          % resting potential (mV)
Eth = -55;            % threshold (mV)
tau = 20;             % membrane time constant (ms)
R = 10;               % membrane resistance (Ohm)
I_inj = 0;            % no constant injected current
sigma_noise = 8.5;  % from part (d)   


% from part (a) i know I-c=1 has 0Hz spiking rate and I-c=2 has 71Hz
% therefore, we can use those as bounds to find I-c that produce spike in
% presence of the noise found in part (d)
I_low=0;
I_high=2;
tolerance=0.01; % i might need to change this later ( difference between high ad low)
min_rate= 0.1;  % i might need to adjust this (min non-zero spiking rate)


nTrials = 5;          %  5 repeated trials
% run multiple trials and then calculate mean firing rate for each I_c
function mean_rate=measure_rate(I_c,sigma_noise,Erest,Eth,tau,R,dt,T,nTrials)
% array for each trial's firing rate
rates= zeros(1,nTrials);
for trial=1:nTrials
    rng(trial); 
     v = Erest * ones(1, round(T/dt)+1);
        spike_count = 0;
        
        % Euler's integration
        for i = 1:length(v)-1
            I_noise = sigma_noise * randn; % random noise
            % updating LIF voltage
            v(i+1) = v(i) + (dt/tau)*( (Erest - v(i)) + R*(I_c + I_noise) );
            
            % threshold crossing check
            if v(i+1) >= Eth
                spike_count = spike_count + 1;
                v(i+1) = Erest; % reset
            end
        end
        
        % Ccalculate the firing rate for this trial (spikes/sec)
        rates(trial) = spike_count / (T/1000); 
    end
    
    % return the mean across all trials
    mean_rate = mean(rates);3
end

%  iteratively narrow down I-low and I-high until they're close
I_candidates = [];
firing_rates_bisect = [];

while (I_high-I_low)>tolerance
    I_mid= (I_low+I_high)/2
    % measure current at the intermediate value 
    current_rate=measure_rate(I_mid,sigma_noise,Erest,Eth,tau,R,dt,T,nTrials)

    I_candidates(end+1) = I_mid;
    firing_rates_bisect(end+1) = current_rate;

    
    fprintf('Testing I_c = %.3f mA -> Avg firing rate = %.2f Hz\n', I_mid, current_rate);
    % now determine whether I_mid produces non-zero firing rate
    % if the rate is is greater than the min-rate we defined, we consider it as
    % spiking
    if current_rate>min_rate
        I_high=I_mid % spike happened so narrow down the upper bound
    else
        I_low=I_mid % no spike so raise up the lower bound
    end
end

% converging threshold is approximately between the updated low and high
new_threshold= (I_low+ I_high)/2;
fprintf('\nEstimated new threshold current (with noise) = %.4f mA\n', new_threshold)


% plotting average firing rate vs. candidate I_c values
figure;
plot(I_candidates, firing_rates_bisect, 'bo-', 'LineWidth', 2);
xlabel('Injected Current I_c (mA)');
ylabel('Average Firing Rate (Hz)');
title('Average Firing Rate vs. Injected Current during Bisection Search');
grid on;
hold on;
% highlight threshold
xline(new_threshold, 'r--', 'LineWidth', 2);
hold off;