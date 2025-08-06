rng(41);
% fundamental parameters
N=10; % 10 trials
r=5; % constant firing rate of poisson process
T=0.1; % total time of simulation
dt=1e-3;


% derived parameters

p= r*dt;
nbins= round(T/dt); % number of bins for each trial
% spike matrix
spk= rand(nbins,N)<p; % bins on rows and trials on columns then compare each value to p
% spk will have 0s(no spike) and 1s(spike) as entries

% firing rate in in each trial (column sum)
num_spk= sum(spk,1);
rate=num_spk/T; % spike rate in Hz

% statistical summary of the spikes
mu=mean(rate);
sigma=std(rate,1);

fprintf(' Mean firing rate m=%.2fHz \n',mu);
fprintf('standard deviation of spikes s=%.2fHz \n',sigma);

% raster plot
figure;  hold on;
for k = 1:N
    t_spk = find(spk(:,k))*dt*1e3; % spike times in ms
    plot(t_spk, k+0*t_spk, 'k|', 'MarkerSize',4,'LineWidth', 2.0)
end
xlabel('time (ms)');  ylabel('trial #');
title(sprintf('Poisson spike trains  (r = %g Hz,  Δt = 1 ms)',r));
xlim([0 T*1e3]);  ylim([0.5 N+0.5]);  box off;





% 
% for trial=1: N
%     spikes= zeros(1,T);
%     for i= 1: T
%         xrand= rand(1);
%         if p> xrand
%             spikes[i]=1;
%         end
%     end
% 
% end

