
%% ----------- problem 3c: Inhomogeneous poisson with exponential refractory period--------------------------------------------------------
% Pool ISIs across M independent 1 s trials for tau_r = 10 ms
dt    = 1e-3;
t     = 0:dt:1;
rfun  = @(tt) 20*sin(4*pi*tt - pi/2) + 30;   % underlying r(t)
tau_r = 10e-3;                              % 10 ms refractory time constant
M     = 100;                                % # of independent trials
all_isis = [];                              % collect ISIs here

rng(0,'twister');
for trial = 1:M
    del_t = Inf;            % time since last spike
    spk   = false(size(t));
    for i = 1:numel(t)
        R_current = rfun(t(i))*(1 - exp(-del_t/tau_r));
        if (R_current*dt) > rand
            spk(i)  = true;
            del_t  = 0;
        else
            del_t = del_t + dt;
        end
    end
    times = find(spk)*dt;         % spike times in seconds
    isis  = diff(times)*1e3;      % ISIs in milliseconds
    all_isis = [all_isis; isis(:)];  %#ok<AGROW>

end

%% --------------------------------------------------------
% Histogram with moderate bin‐width, zooming into [0,50] ms
edges = 0:2:50;   % 2 ms bins
counts = histcounts(all_isis, edges);
prob   = counts / sum(counts);
centers = edges(1:end-1) + diff(edges)/2;

figure('Color','w');
bar(centers, prob, 'FaceColor',[.2 .6 .9]);
hold on;
xline(10,'r--','LineWidth',1.5);
xlim([0 50]);
xlabel('Inter–spike interval \Delta t (ms)');
ylabel('P(\Delta t)');
title(sprintf('Pooled ISI histogram (\\tau_r=%.0f ms)',tau_r*1e3));


















% %% problem 3d: validate 3c with ISI method
% 
% %% ----------From C we have the following------------
% %% problem 3c: Inhomogeneous poisson with exponential refractory period
% dt=1e-3;
% t=0:dt:1;
% nbins=numel(t);
% r= @(t) 20*sin(4*pi*t- pi/2)+30;
% tau_list= [1e-3,10e-3];
% R_all    = zeros(numel(tau_list), nbins); % 1*2 matrix initialized with zeros
% spikes=cell(size(tau_list)); % a 1*2 cell to save spikes for the two taus
% 
% % now let's loop over the refractory constants
% rng(0,'twister'); 
% for k=1:numel(tau_list)
%     tau_r= tau_list(k);
%      del_t    = Inf; % time since last spike
%      Rvec=zeros(1,nbins); % vector of R(t) initialized with zeros for current row of R_all
%      spk=false(1,nbins); % vector of spikes initialized with zeros (logical array)
%      % now let's loop over each time interval for a given refractory constant tau
% 
%      for i=1:nbins
%          % a. calculate r(t) and corresponding R(t)
% 
%          rt= r(t(i)); % current sinusoidal rate
%          R_current= rt*(1-exp(-del_t/tau_r)); % R(t) function of r(t) and current tau
%          Rvec(i)= R_current; % replace with the current R(t)
%          % b. let's draw spike probability with R_current*dt
%          if R_current*dt>rand
%              spk(i)=true;
%              del_t=0; % reset the recovery period
% 
%          else
%              del_t=del_t+dt; % accumulate the recovery period
%          end
% 
%      end
%      R_all(k,:)=Rvec;
%      spikes{k}=find(spk)*dt; % find indices of true enteries in spk. then nbin(i)*dt= t(i)=spike time
% 
% end
% %% -----------------using results from 3c------------------
% %     spikes{k} is a column vector of spike times (in seconds)
% %     for tau_r = tau_list(k).  Here k=2 corresponds to tau_r=10ms.
% M=100;
% all_isis=[];
% % we choose tau=10ms for this case
% for i=1:M
%     spk10=spikes{2};
%     % compute ISIs
% 
%     isis=diff(spk10); % this is in seconds
%     % iss in ms for plot
%     isis_ms= isis*1e3;
%     all_isis=[all_isis;isis];
% end
% 
% edges = 0:2:50;   % 2 ms bins out to 50 ms
% counts = histcounts(all_isis, edges);
% prob   = counts / sum(counts);
% centers = edges(1:end-1) + diff(edges)/2;
% 
% figure('Color','w');
% bar(centers, prob, 'FaceColor',[.2 .6 .9]);
% hold on;
% xline(10,'r--','LineWidth',1.5);
% xlim([0 50]);
% xlabel('Inter–spike interval \Delta t (ms)');
% ylabel('Probability');
% title('Pooled ISI histogram (\tau_r=10\,ms)');
% 
% 
% 
% 
% % % plot of the histogram
% % figure('Color','w');
% % histogram(isis_ms, 'BinWidth',1, 'Normalization','probability');
% % xlabel('Inter–spike interval \Delta t (ms)');
% % ylabel('Probability');
% % title('ISI histogram for tau_r=10ms');
% % 
% % % mark the nominal refractory period
% % hold on;
% % y = ylim;
% % plot(10*[1 1], y, 'r--', 'LineWidth',1.5);
% % text(10, 0.9*y(2), 'tau_r=10ms','Color','black','HorizontalAlignment','left');
% % 
% % % % assume we still have `edges` and `counts` from your histogram code:
% % % centers = edges(1:end-1) + diff(edges)/2;
% % % 
% % % % find all bins with center < tau_r
% % % idx = centers < 10;  
% % % 
% % % figure; 
% % % bar(centers(idx), counts(idx)/sum(counts),'hist');  
% % % xlabel('\Delta t (ms)');  ylabel('Probability');
% % % title('Zoom: ISI < \tau_r = 10 ms');
% % % xlim([0 12]);
% % edges = 0:2:20;                    % 2 ms bins up to 20 ms
% % h = histogram(isis_ms,'BinEdges',edges,...
% %               'Normalization','probability');
% % xlim([0 20])
