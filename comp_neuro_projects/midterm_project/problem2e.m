%% problem2e: reverse‑correlation with “blob” stimuli
% (e) Devise your own STIM design (here: small Gaussian blobs instead of
%     point‐dots) and repeat the RF reconstruction from part d.

clear; clf; rng(0,'twister');

%% 1) define Gabor RF (same as 2a)
sig_x = 5;         % Gaussian sigma along carrier axis (deg)
sig_y = 5;         % Gaussian sigma orthogonal to carrier (deg)
k     = 0.5;       % spatial frequency (cycles/deg)
phi   = 0;         % phase offset
theta = pi/6;      % tilt angle from vertical (rad)

x = -25:25; y = -25:25;
[xx,yy] = meshgrid(x,y);

% rotate coords:
xp =  xx*cos(theta) + yy*sin(theta);
yp = -xx*sin(theta) + yy*cos(theta);

% normalized Gabor envelope
normRF = 1/(2*pi*sig_x*sig_y);
RF = normRF * exp( - (xp.^2)/(2*sig_x^2) - (yp.^2)/(2*sig_y^2) ) ...
         .* cos(k*xp - phi);

%% 2) generate blob‑filtered STIMs
sig_b   = 2;            % blob sigma (deg)
Nwhite  = 100;          % # of +1 blobs per frame
Nblack  = 100;          % # of −1 blobs per frame
Ntrials = 3000;         % total STIM frames
p       = 10;           % keep top p% responses

% build 2D Gaussian blob kernel
w = ceil(3*sig_b);
[xb,yb] = meshgrid(-w:w,-w:w);
B = exp( - (xb.^2+yb.^2)/(2*sig_b^2) ...
         - (yb.^2+yb.^2)/(2*sig_b^2) );
B = B / sum(B(:));

% storage
all_STIM = zeros(numel(y),numel(x),Ntrials);
L        = zeros(Ntrials,1);

%% 3) loop: place random ±1 dots, blur to blobs, renormalize, compute L
for k = 1:Ntrials
  % raw point pattern
  pts = zeros(numel(y),numel(x));
  
  % white centers
  rows = randi(numel(y),[Nwhite,1]);
  cols = randi(numel(x),[Nwhite,1]);
  idx  = sub2ind(size(pts), rows, cols);
  pts(idx) = +1;
  
  % black centers
  rows = randi(numel(y),[Nblack,1]);
  cols = randi(numel(x),[Nblack,1]);
  idx  = sub2ind(size(pts), rows, cols);
  pts(idx) = -1;
  
  % blur to little Gaussian blobs
  s = conv2(pts, B, 'same');
  s = s / max(abs(s(:)));         % renormalize to lie in [−1,1]
  
  all_STIM(:,:,k) = s;
  
  % linear response L_k = ∑_{x,y} RF(x,y) · s(x,y)
  L(k) = sum( RF(:) .* s(:) );
end

%% 4) pick top‑p% by L and reconstruct
thr       = prctile(L,100-p);
keepMask  = L >= thr;
RF_rec    = mean( all_STIM(:,:,keepMask), 3 );

%% 5) display original and reconstructed RF
figure('Color','w','Position',[100 100 800 350]);

subplot(1,2,1)
imagesc(x,y,RF);
axis image xy; colormap(jet); colorbar
title('Original Gabor RF')
xlabel('x (deg)'); ylabel('y (deg)');

subplot(1,2,2)
imagesc(x,y,RF_rec);
axis image xy; colormap(jet); colorbar
title(sprintf('Blob‑STIM RF recon (top %d%%)',p))
xlabel('x (deg)'); ylabel('y (deg)');
