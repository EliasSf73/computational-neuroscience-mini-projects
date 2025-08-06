% Define model parameters
dt=0.1;
Erest= -65;
tau=20;
R=10;
Eth=-55;
t_array=0:dt:1000;
v=zeros(size(t_array));
v(1)=Erest;
spike_times=[];
I_inj=1+1e-13;% try 1+1e-14 to see no spike
% iterate and check 
for i=1:length(t_array)-1
    % check if v(i) crosses the threshold value
    if v(i) >= Eth
        spike_times=[spike_times,t_array(i)];
        v(i+1)=Erest;
    else
        v(i+1)=v(i) + (dt/tau) * (Erest-v(i)+ R*I_inj);
    end


end

% count spikes per second
num_spikes=length(spike_times);
total_time= t_array(end)/1000; % total taken time in seconds
firing_rate=num_spikes/total_time


% % plot
plot(t_array, v);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('voltage trace at I_{inj}=1+1e-14 ')


% Raster plot for the final value of I_inj that produced 22Hz firing rate
figure; 
hold on;
%  plot a vertical line from y=0 to y=1 to represent each spike
for j = 1:length(spike_times)
    line([spike_times(j) spike_times(j)], [0 1], 'Color', 'k', 'LineWidth', 2);
end
xlabel('Time (ms)');
ylabel('Spike');
title('Raster Plot of Spikes (1 sec, I_{inj}=1+1e-13 that yields 1 Hz)');
ylim([0 1.5]); 
hold off;
