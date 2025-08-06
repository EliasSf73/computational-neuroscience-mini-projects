%HW2Prob2a_CoupledLIF.m

%  LIF parameters( units: ms,mV,mA/cm^2)
C=1;
gL=0.1;
tau_m= C/gL;
V_rest=-65;
V_th=-55;
V_reset=-75; % hyperpolarization state
t_ref=2;

T_sim=1000; % simulation for 1second
dt= 0.1;

%  search for Ie that produces 25Hz with bisection method ( for
%  effectiveness)
I_low=0; I_high=2; % initial experimental bounds
stim_onset=0;
target_rate= 25;
tolerance=0.01;
nTrials=5; % multiple trials for averaging

%  a simulation function for one neuron given Ie and returning firing_rate
function rate= lif_rate(Ie,dt, T_sim,gL,V_reset,stim_onset,V_rest,V_th,t_ref,tau_m,nTrials)
numSteps = round(T_sim/dt);
total_spikes=0;

for trial=1:nTrials
    V=V_rest;
    t_last=-Inf;
    spike_count=0;

    for i=1:numSteps   % Loop over each time step
        t_current= i*dt;
        % injecting current at t= stim_onset
        if t_current> stim_onset
            I=Ie;
        else
            I=0;
        end
        % now we have Ie , we'll update the LIF neuron's equation with
        % Euler's method
        % case 1: Neuron is in refractory period
        if (t_current-t_last)<t_ref
            V=V_reset;
        % neuron is not in refractory period: so, we do the update
        else
            V= V+ dt*( (V_rest-V) + (I/gL) )/tau_m;
        end
         % now let's check if the updated neuron crosses threshold
         if V>= V_th && (t_current-t_last)>=t_ref
             spike_count=spike_count+1;
             V=V_reset; % reset V after spike
             t_last=t_current;
         end
    end
    % done with one trial, so we add the spike_count to tota_spikes
    total_spikes=total_spikes+spike_count;
end
% now we have total_spikes for all trials, we'll average them
rate= total_spikes/(nTrials);
rate= rate/((T_sim-stim_onset)/1000); %
end


% Now let's search for Ie that produces 25Hz using bisection method
while (I_high-I_low)> tolerance % an update will continue until they're very close to each other
    I_mid= (I_high+I_low)/2;
    % check rate at this current value
    current_rate=lif_rate(I_mid,dt,T_sim,gL,V_reset,stim_onset,V_rest,V_th,t_ref,tau_m,nTrials);
    fprintf('Testing Ie= %.3f--> average firing rate= %.2f Hz\n', I_mid,current_rate);
    % update will continue until we get the 25 rate and i_high-I_low <=
    % tolerance
    if current_rate< target_rate
        I_low= I_mid;
    else
        I_high=I_mid
    end

end
% now  I_high-I_low <=tolerance: their average would be best solution
I_target= (I_high+I_low)/2;
fprintf('Estimated Ie for 25Hz spiking=%.4f\n',I_target);

% Now we have rate calculating function and Ie that corresponds to 25Hz,
% let's do simulation with Ie and store spike times for raster plot

V=V_rest;
t_vec=0:dt:T_sim;
V_trace=zeros(size(t_vec));
spike_times=[];
t_last=-Inf;
for i=1:length(t_vec)
    % injecting current at t= stim_onset
    t_current=t_vec(i);
    if t_current >= stim_onset
        I=I_target;
    else
        I=0;
    end
    % checking refractory period before updating
    if (t_current-t_last)< t_ref
        V=V_reset;
    else
        V= V+ dt*((V_rest-V) + (I/gL) )/tau_m;
    end
    % save the current V into v vector
    V_trace(i)=V;
    if V>=V_th && (t_current-t_last)>= t_ref
        spike_times(end+1)=t_current;
        V=V_reset;
        t_last=t_current;
    end
end


% plot of the v_traces
figure;
plot(t_vec,V_trace,'LineWidth',2);
xlabel('Time(ms)');
ylabel('Membrane potential(mV)');
title(sprintf('LIF Model V(t) with Ie=%.4f (Target ~25Hz)',I_target));
grid on;


