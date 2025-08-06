%% problem 2b: STM images of N=30,white= +1 and black= -1
clear; clf; rng(0,'twister');
gridSize= 51; % from x=y=-25:25
Nwhite=30; % number of +1 dots
Nblack=30; % number of -1 dots
nsamples=5; % number of independent STMs
figure('Color','W','Position',[100 100 1200 300]);
for k=1:nsamples
    % starting with zero background
    s=zeros(gridSize,gridSize);
    % then let's choose 60 random pixels among 51*51 
    all_ind= randperm(gridSize^2,Nwhite+Nblack);
    white_ind=all_ind(1:Nwhite);
    black_ind= all_ind(Nwhite+1:end); 
    % assign +1 and -1 to white and black respectively
    s(white_ind)=+1;
    s(black_ind)=-1;

    % plot
    subplot(1,nsamples,k);
     imagesc([-25 25],[-25 25],s);
    axis square;
    colormap([0 0 0; 0 1 0; 1 1 1]);  
      % map −1→black, 0→green, +1→white
    caxis([-1 1]);
    set(gca,'YDir','normal');
    title(sprintf('STIM %d',s));
    xlabel('x (deg)'); ylabel('y (deg)');
end



