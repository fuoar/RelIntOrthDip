theta = linspace(0,pi,1e2);
phi=linspace(0,2*pi,1e2);
figure
polarplot(phi,sin(phi).^2)
line([pi 0],[0.5 0.5],'Color','k','LineWidth',3)
ax=gca; ax.ThetaZeroLocation = 'top';
funprep