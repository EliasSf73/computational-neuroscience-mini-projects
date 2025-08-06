c1= 20;
c2= 4*pi;
c3=pi/2;
c4=30;
t= linspace(0,1,1000); % t from 0->1 sec
r=c1*sin(c2*t-c3)+c4;  % r(t)


% plot t vs r(t)
figure('Color','w');
plot(t, r, 'LineWidth', 2);
xlabel('Time t (s)','FontSize',12);
ylabel('Firing rate r(t) (Hz)','FontSize',12);
title('Problem 3a: r(t)=20 sin(4πt - π/2)+30','FontSize',14);
ylim([0 60]);
grid on;