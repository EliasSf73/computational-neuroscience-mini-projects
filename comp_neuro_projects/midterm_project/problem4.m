% Problem4= Problem 2 f: Quantify reconstruction accuracy (ρ) as a function of
% # of STIMs and selection threshold

clear; clf; rng(0,'twister');

%% 1) Define “ground‑truth” Gabor RF from 2a
sig_x = 5;  sig_y = 5;  k = 0.5;  phi = 0;  theta = 30*pi/180;
x = -25:25;  y = -25:25;
[XX,YY] = meshgrid(x,y);
Xp =  XX*cos(theta) + YY*sin(theta);
Yp = -XX*sin(theta) + YY*cos(theta);
% normalized 2‑D Gabor
norm_factor = 1/(2*pi*sig_x*sig_y);
D = norm_factor * exp( - (Xp.^2)/(2*sig_x^2) - (Yp.^2)/(2*sig_y^2) ) ...
            .* cos(k*Xp - phi);

%% 2) Parameters for STIM generation (same as 2b)
Nwhite = 30;    % # white dots
Nblack = 30;    % # black dots
bgval   = 0;    % background

%% 3) Conditions to sweep
trial_counts = [200, 1000, 5000];
pct_list     = [20, 10, 5];   % percentages to keep

%% 4) Pre‑allocate results
rho = zeros(numel(trial_counts), numel(pct_list));

%% 5) Loop over conditions
for ti = 1:numel(trial_counts)
  M = trial_counts(ti);
  
  % 5a) generate M white‐noise STIMs
  STIM = zeros(numel(y), numel(x), M);
  for j = 1:M
    pts = zeros(numel(y), numel(x));
    % place +1 dots
    idx = randperm(numel(pts), Nwhite);
    pts(idx) = +1;
    % place -1 dots (avoid collisions)
    empty = find(pts==0);
    idx2  = empty(randperm(numel(empty), Nblack));
    pts(idx2) = -1;
    STIM(:,:,j) = pts;
  end
  
  % 5b) compute linear responses L_j = sum(D.*STIM)
  Dj = D(:);
  L  = reshape( Dj' * reshape(STIM,[],M), M, 1 );
  
  for pi = 1:numel(pct_list)
    p = pct_list(pi);
    thresh = prctile(L, 100 - p);
    keep  = L > thresh;
    
    % 5c) reconstruct RF: average of selected STIMs
    Drec = mean( STIM(:,:,keep), 3 );
    
    % 5d) compute Pearson correlation ρ between D and Drec
    Dr = D(:);       Drm = mean(Dr);
    Rec = Drec(:);   Recm = mean(Rec);
    rho_val = ( (Dr-Drm)'*(Rec-Recm) ) ...
              / ( norm(Dr-Drm)*norm(Rec-Recm) );
    rho(ti,pi) = rho_val;
  end
end

%% 6) Display results
fprintf('Reconstruction ρ for different #STIMs (rows) and %%keep (cols):\n');
disp(array2table(rho, 'VariableNames', ...
    strcat(string(pct_list),'%'), ...
    'RowNames', strcat(string(trial_counts),' trials')) );

%% 7) Plot heatmap
figure('Color','w');
imagesc(pct_list, trial_counts, rho);
set(gca,'YDir','normal','XDir','reverse');
xlabel('Selection threshold (%)');
ylabel('Number of STIMs');
title('rho: corr(D,hatD) for varying trials & %keep');
colorbar;
