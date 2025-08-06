%% problem 3c: Inhomogeneous poisson with exponential refractory period
dt=1e-3;
t=0:dt:1;
nbins=numel(t);
r= @(t) 20*sin(4*pi*t- pi/2)+30;
tau_list= [1e-3,10e-3];
R_all    = zeros(numel(tau_list), nbins); % 1*2 matrix initialized with zeros
spikes=cell(size(tau_list)); % a 1*2 cell to save spikes for the two taus

% now let's loop over the refractory constants
rng(0,'twister'); 
for k=1:numel(tau_list)
    tau_r= tau_list(k);
     del_t    = Inf; % time since last spike
     Rvec=zeros(1,nbins); % vector of R(t) initialized with zeros for current row of R_all
     spk=false(1,nbins); % vector of spikes initialized with zeros (logical array)
     % now let's loop over each time interval for a given refractory constant tau

     for i=1:nbins
         % a. calculate r(t) and corresponding R(t)

         rt= r(t(i)); % current sinusoidal rate
         R_current= rt*(1-exp(-del_t/tau_r)); % R(t) function of r(t) and current tau
         Rvec(i)= R_current; % replace with the current R(t)
         % b. let's draw spike probability with R_current*dt
         if R_current*dt>rand
             spk(i)=true;
             del_t=0; % reset the recovery period
         
         else
             del_t=del_t+dt; % accumulate the recovery period
         end

     end
     R_all(k,:)=Rvec;
     spikes{k}=find(spk)*dt; % find indices of true enteries in spk. then nbin(i)*dt= t(i)=spike time

end

% plot
figure('Color','w');
plot(t, r(t),      'k-', 'LineWidth',1.5), hold on;
plot(t, R_all(1,:), 'b-');
plot(t, R_all(2,:), 'r-');
legend('r(t)','R_{1\rm ms}(t)','R_{10\rm ms}(t)','Location','NorthEast');
xlabel('Time (s)'), ylabel('Rate (Hz)');
title('Problem 3c: Refractory‐modulated rate');





