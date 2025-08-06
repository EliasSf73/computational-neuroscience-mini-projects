%% ---------------------------------------------------------------
%  check_minTime.m  –  verifies §1 d  “400‑spike” rule–of–thumb
%  
% ---------------------------------------------------------------
clear;  clc;  close all;

%% DESIGN TARGETS --------------------------------------------------
eps_r  = 0.05;           % 5 % precision for rate (relative SD)
eps_F  = 0.10;           % 10 % SD for Fano‑factor
Ntr    = 30;             % trials per condition
rates  = 5:5:50;         % Hz  (edit freely)

M      = 1e4;            % Monte‑Carlo repetitions
rng(1)                   % reproducible

%% MINIMUM #SPIKES → TIME PER TRIAL -------------------------------
nMin_rate = ceil(1/eps_r^2);                 % from Var[\hat r]/r^2
nMin_Fano = ceil(2/((Ntr-1)*eps_F^2));       % from Var[\hat F]
nMin      = max(nMin_rate, nMin_Fano);       % overall
Tmin      = nMin ./ rates;                  % s  per trial

fprintf('Rule‑of‑thumb: need %d spikes  (%.0f trials)\n',nMin,Ntr)

%% MONTE‑CARLO VERIFICATION ---------------------------------------
relSD_rate = zeros(size(rates));
SD_Fano    = zeros(size(rates));

% choose built‑in poissrnd if available, otherwise use fallback
hasPRND = exist('poissrnd','file')==2;
if ~hasPRND,  prnd = @poissrnd_local;  else,  prnd = @poissrnd;  end

for k = 1:numel(rates)
    r   = rates(k);                % Hz
    T   = Tmin(k);                 % s
    lam = r*T;                     % Poisson mean / trial

    counts   = prnd(lam,M,Ntr);    % [M × Ntr] spike counts
    mu_hat   = mean(counts,2)/T;                 % Hz
    var_hat  = var(counts,0,2)/(T^2);            % Hz²
    Fano_hat = var_hat./mu_hat;

    relSD_rate(k) = std(mu_hat)/mean(mu_hat);    % SD/mean
    SD_Fano(k)    = std(Fano_hat);               % SD(F̂)
end

%% NUMERICAL TABLE -------------------------------------------------
tbl = table(rates(:),Tmin(:),relSD_rate(:),SD_Fano(:), ...
     'VariableNames',{'Rate_Hz','Tmin_s','RelSD_rate','SD_Fano'});
disp(tbl)

%% FIGURES ---------------------------------------------------------
figure('Color','w');

subplot(1,2,1)
plot(rates,relSD_rate,'ko-','MarkerFaceColor','k'); hold on
yline(eps_r,'r--','LineWidth',1.4);
xlabel('rate  r  (Hz)');
ylabel('relative SD(\hat{\mu})','Interpreter','latex');
title('Rate–estimate precision');
grid on

subplot(1,2,2)
plot(rates,SD_Fano,'ko-','MarkerFaceColor','k'); hold on
yline(eps_F,'r--','LineWidth',1.4);
xlabel('rate  r  (Hz)');
ylabel('SD(\hat{F})','Interpreter','latex');
title('Fano‑factor precision');
grid on

%% --------- local Poisson generator (fallback) -------------------
function X = poissrnd_local(lambda,m,n)
% Very small λ  →  Knuth algorithm
% Large  λ      →  normal approximation
if lambda < 30
    L = exp(-lambda);
    X = zeros(m,n);
    for j = 1:n
        for i = 1:m
            k = 0;  p = 1;
            while p > L
                k = k + 1;
                p = p * rand;
            end
            X(i,j) = k-1;
        end
    end
else
    X = round( lambda + sqrt(lambda).*randn(m,n) );
    X(X<0) = 0;
end
end
