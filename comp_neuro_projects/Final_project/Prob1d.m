%% Problem 1d: STDP Learning Window
% Compute and plot ΔEPSC (%) vs. Δt = t_post – t_pre over [-100,100] ms
% using exponential STDP with τ₊ = τ₋ = 10 ms and amplitudes ±100%.

clear; close all;

%% 1) STDP parameters
A_plus  = 100;    % maximum potentiation (%)
A_minus = 100;    % maximum depression (%)
tau_plus  = 10;   % LTP time constant (ms)
tau_minus = 10;   % LTD time constant (ms)

%% 2) Timing vector
dt       = 1;                % resolution (ms)
dT       = -100:dt:100;      % Δt = t_post - t_pre
nPoints  = numel(dT);

%% 3) Compute ΔEPSC (%) for each Δt
deltaEPSC = zeros(size(dT));
for i = 1:nPoints
    if dT(i) > 0
        % post after pre → potentiation
        deltaEPSC(i) = A_plus * exp(-dT(i)/tau_plus);
    elseif dT(i) < 0
        % pre after post → depression
        deltaEPSC(i) = -A_minus * exp(dT(i)/tau_minus);
    else
        % Δt = 0, peak LTP
        deltaEPSC(i) =  A_plus;
    end
end

%% 4) Plot the STDP window
figure('Name','STDP Learning Window','NumberTitle','off');
plot(dT, deltaEPSC, 'LineWidth',2);
xlabel('\Delta t = t_{post} - t_{pre} (ms)', 'FontSize',12);
ylabel('\Delta EPSC (%)', 'FontSize',12);
title('STDP Window: Exponential LTP / LTD (\tau = 10 ms)', 'FontSize',14);
grid on;


