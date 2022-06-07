figure
n1 = 1.5; n2 = 1;
NA = 0.5;
theta2lim = asin(NA./n2);
thetas2lim = linspace(0,theta2lim,1e2);
thetas2 = linspace(0,pi/2,1e2);
fn21 = @(x) rad2deg(asin(n2./n1.*sin(x)));
plot(rad2deg(thetas2), fn21(thetas2),'DisplayName','Snell-Descartes')
hold on
for NA=[0.2:0.2:1]
%NA = .65;
theta2lim = asin(NA./n2);
thetas2lim = linspace(0,theta2lim,1e2);
thetas2 = linspace(0,pi/2,1e2);
fn21 = @(x) rad2deg(asin(n2./n1.*sin(x)));
%plot(rad2deg(thetas2), fn1(thetas2))
hold on
plot([1 1].*max(fn21(thetas2lim)),[0 40],'DisplayName',sprintf('NA = %.1f',NA))
end
plot([1 1].*rad2deg(asin(n2./n1)),[0 40],'r--','DisplayName','TIR')
legend('Location','best')
title(sprintf('Dipole in medium 1 n_1 = %.1f collected in medium 2 n_2 = %.1f', n1, n2))
xlabel('Refracted angles \theta_2 (°)'); ylabel('Corresponding incidence angle \theta_1 (°)')





%%
fn12 = @(x) rad2deg(asin(n1./n2.*sin(x)));
thetas1 = linspace(0,pi/2,1e2);

figure
plot(rad2deg(thetas1), fn12(thetas1),'DisplayName','1->2')
hold on
%plot(rad2deg(thetas1), fn21(thetas1),'DisplayName','2->1')
legend
for NA=[0.25:0.25:1 0.65]
%NA = .65;
theta2lim = asin(NA./n2);
thetas2lim = linspace(0,theta2lim,1e2);
thetas2 = linspace(0,pi/2,1e2);
fn21 = @(x) rad2deg(asin(n2./n1.*sin(x)));
%plot(rad2deg(thetas2), fn1(thetas2))
hold on
plot([1 1 0].*max(fn21(thetas2lim)),[0 rad2deg(theta2lim) rad2deg(theta2lim)],'DisplayName',sprintf('NA = %.2f',NA))
%plot([0 max(fn21(thetas2lim))],[1 1].*rad2deg(theta2lim),'DisplayName',sprintf('NA = %.1f',NA))

end

theta1lim = rad2deg(asin(n2./n1));
plot([1 1].*theta1lim, [0 90],'b-.','DisplayName','TIR')
xlabel('Incidence angle \theta_{in}'); ylabel('Refracted angle \theta_{r}');

title(sprintf('Medium 1 n_1 = %.2f, medium 2 n_2 = %.2f', n1, n2))
