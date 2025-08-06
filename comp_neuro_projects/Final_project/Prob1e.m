%% Problem 1e: Estimate ΔEPSC After a Single Bursting Input via STDP
% This standalone script:
% 1) Generates the bursting presynaptic spike train (Part 1b parameters)
% 2) Simulates the postsynaptic LIF response (using best_w from Part 1c)
% 3) Applies the exponential STDP rule (τ₊=τ₋=10 ms, A₊=A₋=100%) to compute
%    the net percent change in synaptic efficacy ΔEPSC_total
% 4) Updates the EPSC kernel amplitude and plots “before” vs. “after”

clear; close all;

%% 1) Simulation parameters
dt     = 1;              % ms
T      = 500;            % ms
t      = 0:dt:T;         
N      = numel(t);
I_drive = 10;            % tonic drive for bursting model (mA)

%% 2) Bursting presynaptic spike train (Izhikevich parameters)
a = 0.02; b = 0.15; c = -50; d = 2;
V = -65*ones(1,N);
u = b*V;
spikes_bu = [];
for i = 2:N
    dV = 0.04*V(i-1)^2 + 5*V(i-1) + 140 - u(i-1) + I_drive;
    du = a*(b*V(i-1) - u(i-1));
    V(i) = V(i-1) + dt*dV;
    u(i) = u(i-1) + dt*du;
    if V(i) >= 30
        spikes_bu(end+1) = t(i);  %#ok<*AGROW>
        V(i) = c;
        u(i) = u(i) + d;
    end
end

%% 3) Postsynaptic LIF response using optimal weight from 1c
best_w  = 1.09;           % chosen weight from part 1c
tau_syn = 10;             % ms
[spikes_post, ~] = run_LIF_with_syn(best_w, tau_syn, spikes_bu);

%% 4) STDP parameters and net ΔEPSC computation
A_plus   = 100;           % % potentiation at Δt=0
A_minus  = 100;           % % depression at Δt=0
tau_plus  = 10;           % LTP time constant (ms)
tau_minus = 10;           % LTD time constant (ms)

deltaEPSC_total = 0;
for i = 1:numel(spikes_bu)
    del_t = spikes_post - spikes_bu(i);       % vector of post–pre intervals
    % potentiation: Δt > 0
    pos = del_t > 0;
    deltaEPSC_total = deltaEPSC_total + sum( A_plus  * exp(-del_t(pos)/tau_plus) );
    % depression: Δt < 0
    neg = del_t < 0;
    deltaEPSC_total = deltaEPSC_total - sum( A_minus * exp( del_t(neg)/tau_minus) );
end

fprintf('Net ΔEPSC after one burst: %.1f%%\n', deltaEPSC_total);

%% 5) Update EPSC kernel amplitude and plot
w_before = best_w;
w_after  = w_before * (1 + deltaEPSC_total/100);

t_syn = 0:1:200;                   % ms for EPSC kernel
I_pre  = w_before * exp(-t_syn/tau_syn);
I_post = w_after  * exp(-t_syn/tau_syn);

figure('Name','EPSC Before vs After Burst','NumberTitle','off');
plot(t_syn, I_pre, '--k','LineWidth',1.5); hold on;
plot(t_syn, I_post,'-r','LineWidth',1.5);
hold off;
xlabel('Time (ms)');
ylabel('EPSC (arb. units)');
title(sprintf('EPSC Before/After Burst (ΔEPSC = %.1f%%)', deltaEPSC_total));
legend('Before','After','Location','Best');
grid on;

%% --- Helper function: LIF + exponential synapse ---
function [spike_times_post, V] = run_LIF_with_syn(w, tau_syn, pre_spikes)
  dt = 1; T = 500; t = 0:dt:T; N = numel(t);
  V = -65 * ones(1,N);
  I_syn = zeros(1,N);
  spike_times_post = [];
  idx_pre = 1;

  for i = 2:N
    % decay synaptic current
    I_syn(i) = I_syn(i-1) - (dt/tau_syn)*I_syn(i-1);
    % inject at presyn spike times
    if idx_pre <= numel(pre_spikes) && t(i)==pre_spikes(idx_pre)
      I_syn(i) = I_syn(i) + w;
      idx_pre = idx_pre + 1;
    end
    % LIF update
    dV = -(V(i-1)+65)/20 + I_syn(i)/1;
    V(i) = V(i-1) + dt*dV;
    if V(i) >= -40
      spike_times_post(end+1) = t(i);
      V(i) = -65;
    end
  end
end
