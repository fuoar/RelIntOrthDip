clear all
close all
clc
clrmp = @(x) brewermap(x,"PuOr");

% Code to modelise the emission pattern of a dipole
% in free space depending on the observation direction


%% Parameters

% dipole orientation
theta_dip = pi/2;
phi_dip = pi/2;

n1 = 1.5;	% PS or SiO2
n2 = 1;		% air

f = 2e-3; % in m
%% Create variables

u_dip = [sin(theta_dip)*cos(phi_dip); ...
	sin(theta_dip)*sin(phi_dip); ...
	cos(theta_dip)]; % dipole unit vector 

% observation direction
thetas_obs = linspace(0,pi,1e2);
phis_obs = linspace(0,2*pi,1e2);

% 
% % transmission coefficients at interfaec
% t_s = (2*n2*cos(theta_obs))./(n1*cos(theta_obs)+n_2*cos(theta_obs2));
% t_p = (2*n1*cos(theta_obs))./(n1*cos(theta_obs2)+n_2*cos(theta_obs));

% electric field from definitions
Edef_tot = zeros(length(thetas_obs),length(phis_obs));	% norme du champ electrique
Edef_s = zeros(length(thetas_obs),length(phis_obs));	% norme du champ selon s
Edef_p = zeros(length(thetas_obs),length(phis_obs));	% norme du champ selon p

% expression Lethiec
Emat_tot = zeros(length(thetas_obs),length(phis_obs));
Emat_s = zeros(length(thetas_obs),length(phis_obs));	% norme du champ selon s
Emat_p = zeros(length(thetas_obs),length(phis_obs));	% norme du champ selon p
Emat_3tot = zeros(length(thetas_obs),length(phis_obs),3);

Pol = zeros(length(thetas_obs),length(phis_obs),1e2);

%% loop over all observation directions
for i=1:length(thetas_obs)
	for j=1:length(phis_obs)
		theta_obs = thetas_obs(i);
		phi_obs = phis_obs(j);


theta_obs2 = asin(n2./n1.*sin(theta_obs));
		% unit vector in observation direction
		u_obs = [sin(theta_obs)*cos(phi_obs); ...
			sin(theta_obs)*sin(phi_obs); ...
			cos(theta_obs)];

		%%% basic definition of electric field as u_obs x u_dip x u_obs
		Eraw = cross(u_obs, cross(u_dip,u_obs));

		% define polarisation unit vectors
		u_s = [-sin(phi_obs); ...
			cos(phi_obs); ...
			0];
		u_p = [cos(theta_obs)*cos(phi_obs);...
			cos(theta_obs)*sin(phi_obs);...
			-sin(theta_obs)];

		%%% define electric field components as projections on unit S and P
		E_s = dot(Eraw, u_s);
		E_p = dot(Eraw, u_p);

		%%% expressions from Lethiec
		E_s_mat = sin(theta_dip)*sin(phi_obs-phi_dip);
		E_p_mat = - sin(theta_obs)*cos(theta_dip) + ...
 			cos(theta_obs)*sin(theta_dip)*cos(phi_dip-phi_obs);

		%%%% from def
		Edef_tot(i,j) = norm(Eraw); % get value
		Edef_s(i,j) = norm(E_s*u_s);
		Edef_p(i,j) = norm(E_p*u_p);

		%%%% from Lethiec
		E_tot_fs = E_s_mat*u_s+E_p_mat*u_p;

		Emat_tot(i,j) = norm(E_tot_fs);
		Emat_3tot(i,j,:) = E_tot_fs;
		Emat_s(i,j) = norm(E_s_mat*u_s);
		Emat_p(i,j) = norm(E_p_mat*u_p);

%%%%%
% after objective
v_s = u_s;

v_p = [cos(phi_obs);...
	sin(phi_obs);...
	0];

rho = f*sin(theta_obs);

Eobj = sqrt(n2./cos(theta_obs2)).*(dot(E_tot_fs,v_s)*v_s+dot(E_tot_fs,v_p)*v_p);
%size(Eobj)
	
		%%% projection on polariser axis

	%	figure(4)
		%if j<i
		%Pol = [];
		alpha = linspace(0,2*pi,1e2);
		for kk=1:length(alpha)%alpha=linspace(0,2*pi,1e2)
		u_alpha = [ cos(alpha(kk));...
			sin(alpha(kk));...
			0];

		%%% from def
		Pol(i,j,kk) = dot(Eobj,u_alpha).^2;
		end


%polarplot(linspace(0,2*pi,1e2),Pol)
%hold on
		%end

	
	
	end
	DOP(i) = getDOP(Emat_tot(i,:));
end



% figure
%% Store simulation results

Simu.Date = sprintf('%s',datetime);
%%% Store simulation details
Simu.thetaphi_dip = [theta_dip; phi_dip];
Simu.thetaphi_obs = [thetas_obs; phis_obs];

Simu.Edef = Edef_tot;
Simu.Edef_s = Edef_s;
Simu.Edef_p = Edef_p;


Simu.Emat = Emat_tot;
Simu.Emat_s = Emat_s;
Simu.Emat_p = Emat_p;


Simu.Pol = Pol;


clearvars -except Simu clrmp

%% Show calculated results


[x,y,z]=PlotResults(Simu,'Lethiec',clrmp);

Simu.dip.x = x; 
Simu.dip.y = y;
Simu.dip.z = z;

%% Function
clearvars -except Simu


function [x,y,z] = PlotResults(Simu,mode,clrmp)

switch mode
	case 'Lethiec'
thetas_obs = Simu.thetaphi_obs(1,:);
phis_obs = Simu.thetaphi_obs(2,:);

theta_dip = Simu.thetaphi_dip(1);
phi_dip = Simu.thetaphi_dip(2);

Edef_tot = Simu.Emat;

Edef_s = Simu.Emat_s;
Edef_p = Simu.Emat_p;

	case 'fromdef'
thetas_obs = Simu.thetaphi_obs(1,:);
phis_obs = Simu.thetaphi_obs(2,:);

theta_dip = Simu.thetaphi_dip(1);
phi_dip = Simu.thetaphi_dip(2);

Edef_tot = Simu.Edef;

Edef_s = Simu.Edef_s;
Edef_p = Simu.Edef_p;

end
% 
% figure
% patternCustom(Edef_tot',rad2deg(thetas_obs), rad2deg(phis_obs))

% figure
% nexttile
% imagesc(rad2deg(Simu.thetaphi_obs(1,:)),rad2deg(Simu.thetaphi_obs(2,:)),Simu.Edef')
% colorbar; caxis([0 1])
% nexttile
% imagesc(rad2deg(Simu.thetaphi_obs(1,:)),rad2deg(Simu.thetaphi_obs(2,:)),Simu.Edef_s')
% colorbar; caxis([0 1])
% nexttile
% imagesc(rad2deg(Simu.thetaphi_obs(1,:)),rad2deg(Simu.thetaphi_obs(2,:)),Simu.Edef_p')
% colorbar; caxis([0 1])

figure('Position',[712.2000 49 1.3520e+03 1.0832e+03], 'Color','w');
tiledlayout('flow','TileSpacing','compact','Padding','compact');

% 
% azel = phitheta2azel([rad2deg(phis_obs);rad2deg(thetas_obs)],false);
% 
nexttile([2 2])
[x,y,z] = spherical2cart(Edef_tot',thetas_obs,phis_obs);
s=surf(x,y,z,Edef_tot);
s(1).EdgeColor='none'; colorbar; colormap(clrmp([]))
pbaspect([1 1 1])
% 
% nexttile
% p=pcolor(rad2deg(thetas_obs)-rad2deg(theta_dip),rad2deg(phis_obs)-rad2deg(phi_dip),Edef_tot');
% p.EdgeColor = 'none';
% colorbar; caxis([0 1]);
% xlabel('\Theta_{obs} - \Theta_{dip}','Interpreter','tex');
% ylabel('\Phi_{obs} - \Phi_{dip}','Interpreter','tex');
% 
% title('Norm of E vs \Theta_{obs} and \Phi_{obs}','interpreter','tex')
% subtitle(sprintf('\\Theta_{dip} = %d and \\Phi_{dip} = %d', rad2deg(theta_dip), rad2deg(phi_dip)));%,'HorizontalAlignment','left');
% 
% ax=gca; ax.TitleHorizontalAlignment = 'left'; 
% %colormap(brewermap([],"PuOr"))
% 
% 
% nexttile
% ss=surfc(rad2deg(thetas_obs), rad2deg(phis_obs), Edef_tot');
% ss(1).EdgeColor='none';
% xlabel('\Theta_{obs}','Interpreter','tex');
% ylabel('\Phi_{obs}','Interpreter','tex');
% zlabel('|E|')



nexttile
for i=1:5:size(Edef_tot,2)
polarplot(thetas_obs-theta_dip,Edef_tot(:,i));%, 'Color',cg(i,:))
if i==1; hold on; end
end
thetalim([0 180]-rad2deg(theta_dip)); title('|E| = f(\Theta_{obs}-\Theta_{dip}) varying \Phi_{obs}');
ax=gca;
ll = length(ax.Children); cg = clrmp(ll);%cgrad3([140 81 10]./255, [0.8 0.8 0.8], [1 102 94]./255, ll);
ax.ColorOrder = cg;

nexttile
for i=1:5:size(Edef_tot,1)
polarplot(phis_obs,Edef_tot(i,:));
if i==1; hold on; end
end
 title('|E| = f(\Phi_{obs}) varying \Theta_{obs}')
ax=gca;
ll = length(ax.Children); cg = clrmp(ll);%cgrad3([140 81 10]./255, [0.8 0.8 0.8], [1 102 94]./255, ll);
ax.ColorOrder = cg;


nexttile
for i=1:5:size(Edef_s,2)
polarplot(thetas_obs,Edef_s(:,i));%, 'Color',cg(i,:))
if i==1; hold on; end
end
thetalim([0 180]); title('|E_{s\Phi}| = f(\Theta_{obs})');
ax=gca;
ll = length(ax.Children); cg = clrmp(ll);%cgrad3([140 81 10]./255, [0.8 0.8 0.8], [1 102 94]./255, ll);
ax.ColorOrder = cg;

nexttile
for i=1:5:size(Edef_s,1)
polarplot(phis_obs,Edef_s(i,:));
if i==1; hold on; end
end
 title('|E_{s\Phi}| = f(\Phi_{obs})')
ax=gca;
ll = length(ax.Children); cg = clrmp(ll);%cgrad3([140 81 10]./255, [0.8 0.8 0.8], [1 102 94]./255, ll);
ax.ColorOrder = cg;


nexttile
funprep

nexttile
for i=1:5:size(Edef_p,2)
polarplot(thetas_obs,Edef_p(:,i));%, 'Color',cg(i,:))
if i==1; hold on; end
end
thetalim([0 180]); title('|E_{p\Theta}| = f(\Theta_{obs})');
ax=gca;
ll = length(ax.Children); cg = clrmp(ll);%cgrad3([140 81 10]./255, [0.8 0.8 0.8], [1 102 94]./255, ll);
ax.ColorOrder = cg;

nexttile
for i=1:5:size(Edef_p,1)
polarplot(phis_obs,Edef_p(i,:));
if i==1; hold on; end
end
 title('|E_{p, u_\Theta}| = f(\Phi_{obs})')
ax=gca;
ll = length(ax.Children); cg = clrmp(ll);%cgrad3([140 81 10]./255, [0.8 0.8 0.8], [1 102 94]./255, ll);
ax.ColorOrder = cg;

sgtitle(sprintf('Simulation of E from dipole with \\Theta = %.0f° and \\Phi = %.0f°: %s',round(rad2deg(theta_dip)), round(rad2deg(phi_dip)),mode))
end

function DOP = getDOP(Evsphi)
	res = (max(Evsphi)-min(Evsphi))./(max(Evsphi)+min(Evsphi));
DOP = res;
end


function [x, y, z] = spherical2cart(r, theta,phi)
x = zeros(size(r,1),size(r,2));
y = zeros(size(r,1),size(r,2));
z = zeros(size(r,1),size(r,2));

for i=1:length(x)
	for j=1:length(y)
	x(i,j) = r(j,i)*sin(theta(i))*sin(phi(j));
	y(i,j) = r(j,i)*sin(theta(i))*cos(phi(j));
	z(i,j) = r(j,i)*cos(theta(i));
	end
end
end