clear all
%close all
clc
clrmp = @(x) brewermap(x,"PuOr");

% Code to modelise the emission pattern of a dipole
% in free space depending on its orientation


%% Parameters

% dipole orientation
theta_dip = pi/2;
phi_dip = 0;


% vary theta dipole
theta_dips = linspace(0,pi/2,1e1);

DOPDOP = [];
figure
for theta_dip=theta_dips

NA = 0.6;

n1 = 1.5;	% PS or SiO2
n2 = 1;		% air

f = 1e-2; % in m 
%% Create variables

%%% by definition we are looking from z+
thistheta_obs = 0;

%%% from NA wa can calculate theta1max
theta2lim = asin(NA./n2);
theta1lim = asin(n2/n1.*sin(theta2lim));


u_dip = [sin(theta_dip)*cos(phi_dip); ...
	sin(theta_dip)*sin(phi_dip); ...
	cos(theta_dip)]; % dipole unit vector

% observation direction
thetas_obs = linspace(0,pi,5e1);
phis_obs = linspace(0,2*pi,1e2);



% electric field from definitions
Edef_tot = zeros(length(thetas_obs),length(phis_obs));	% norme du champ electrique
Edef_s = zeros(length(thetas_obs),length(phis_obs));	% norme du champ selon s
Edef_p = zeros(length(thetas_obs),length(phis_obs));	% norme du champ selon p

% expression Lethiec
Emat_tot = zeros(length(thetas_obs),length(phis_obs));
Emat_s = zeros(length(thetas_obs),length(phis_obs));	% norme du champ selon s
Emat_p = zeros(length(thetas_obs),length(phis_obs));	% norme du champ selon p
Emat_3tot = zeros(length(thetas_obs),length(phis_obs),3);

Pol = zeros(length(thetas_obs(thetas_obs<=theta1lim)),length(phis_obs),1e2);
%DOP = zeros(length(thetas_obs(thetas_obs<=theta1lim)),length(phis_obs));

%% loop over all observation directions to get emission pattern
for i=1:length(thetas_obs)
	for j=1:length(phis_obs)
		%%%%% Get angles and observation unit vector
		theta_obs = thetas_obs(i);
		phi_obs = phis_obs(j);
		theta_obs2 = asin(n2./n1.*sin(theta_obs)); % Snell to get theta_obs2

		% unit vector in observation direction
		u_obs = [sin(theta_obs)*cos(phi_obs); ...
			sin(theta_obs)*sin(phi_obs); ...
			cos(theta_obs)];


		%%%%%
		%%%%% Calculate field as: u_obs x u_dip x u_obs
		Eraw = cross(u_obs, cross(u_dip,u_obs));

		%%%%% Define polarisation unit vectors in plane perp. to u_obs
		u_s = [-sin(phi_obs); ...
			cos(phi_obs); ...
			0];
		u_p = [cos(theta_obs)*cos(phi_obs);...
			cos(theta_obs)*sin(phi_obs);...
			-sin(theta_obs)];

		%%%%% Get projections of field on S and P
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
		%%%%% Interface from n1 to n2
		%%%%%
		% transmission coefficients at interfaec
		t_s = (2*n2*cos(theta_obs))./(n1*cos(theta_obs)+n2*cos(theta_obs2));
		t_p = (2*n1*cos(theta_obs))./(n1*cos(theta_obs2)+n2*cos(theta_obs));

		%%%% Updated field expression with Fresnel coeffients
		F = (E_s*t_s*u_s + E_p*t_p*u_p)./f;
		Edef_after_interface = F.*n2./n1.*cos(theta_obs2)./cos(theta_obs);

		%%%% Updated polarisation unit vectors after objective
		v_s = u_s;
		v_p = [cos(phi_obs);...
			sin(phi_obs);...
			0];

		%%%% Updated field after objective
		Edef_after_obj = E_s*v_s + E_p*v_p;
		rho = f*sin(theta_obs2);
		Eobj = sqrt(n2./cos(theta_obs2)).*Edef_after_obj;

		%%%%%
		%%%%% Projection on polariser
		alpha = linspace(0,2*pi,1.5e2);
		if theta_obs<theta1lim %%% for all angles within NA
			%%% projection on polariser axis
			for kk=1:length(alpha)%alpha=linspace(0,2*pi,1e2)
				u_alpha = [ cos(alpha(kk));...
					sin(alpha(kk));...
					0];

				%%% from def
				rres = dot(Eobj,u_alpha).^2.*rho;
				

				Pol(i,j,kk) = rres;
			end
%getDOP(Pol(i,j,:)))
 				%Imax = max(squeeze(Pol(i,j,:)))
 				%Imin = min(squeeze(Pol(i,j,:)))
				%DOPDOP = [DOPDOP (Imax-Imin)./(Imax+Imin)];
			%DOP(i,j) = getDOP(squeeze(Pol(i,j,:)));
		end%

	end
	
end

SumPolars = squeeze(sum(squeeze(sum(Pol(:,:,:),1)),1));
DOPDOP = [DOPDOP getDOP(SumPolars)];

% figure
%% Store simulation results

Simu.Date = sprintf('%s',datetime);
%%% Store simulation details
Simu.thetaphi_dip = [theta_dip; phi_dip];
Simu.theta_obs = thetas_obs;
Simu.phi_obs = phis_obs;

Simu.Edef = Edef_tot;
Simu.Edef_s = Edef_s;
Simu.Edef_p = Edef_p;


Simu.Emat = Emat_tot;
Simu.Emat_s = Emat_s;
Simu.Emat_p = Emat_p;


Simu.Pol = Pol;
Simu.DOP = DOPDOP;



%% Show calculated results


%[x,y,z]=PlotResults(Simu,'fromdef',clrmp);

%Simu.dip.x = x;
%Simu.dip.y = y;
%Simu.dip.z = z;

%%

%figure('Color','w')
subplot(121)
yy=squeeze(sum(squeeze(sum(Simu.Pol(:,:,:),1)),1));
polarplot(alpha,yy)%./max(max(yy)))
hold on
title(sprintf('\\Theta_{dip} - \\Theta_{obs} = %.1f°', rad2deg(theta_dip-thistheta_obs)))

subplot(122)
polarplot(alpha,yy./max(max(yy)))
hold on
title(sprintf('\\Theta_{dip} - \\Theta_{obs} = %.1f°', rad2deg(theta_dip-thistheta_obs)))

end

ax=gca; ll = length(ax.Children); ax.ColorOrder = clrmp(ll);



figure
nexttile
plot(rad2deg(theta_dips),DOPDOP)
%% Function

%clearvars -except Simu clrmp


function [x,y,z] = PlotResults(Simu,mode,clrmp)

switch mode
	case 'Lethiec'
		thetas_obs = Simu.theta_obs;
		phis_obs = Simu.phi_obs;
		theta_dip = Simu.thetaphi_dip(1);
		phi_dip = Simu.thetaphi_dip(2);
		Edef_tot = Simu.Emat;
		Edef_s = Simu.Emat_s;
		Edef_p = Simu.Emat_p;
	case 'fromdef'
		thetas_obs = Simu.theta_obs;
		phis_obs = Simu.phi_obs;
		theta_dip = Simu.thetaphi_dip(1);
		phi_dip = Simu.thetaphi_dip(2);
		Edef_tot = Simu.Edef;
		Edef_s = Simu.Edef_s;
		Edef_p = Simu.Edef_p;

end
figure('Position',[712.2000 49 1.3520e+03 1.0832e+03], 'Color','w');
tiledlayout('flow','TileSpacing','compact','Padding','compact');


nexttile([2 2])
[x,y,z] = spherical2cart(Edef_tot.^2,thetas_obs,phis_obs);
s=surf(x,y,z,Edef_tot);
s(1).EdgeColor='none';
colorbar; colormap(clrmp([]))
pbaspect([1 1 1]);xlim([-1 1]); ylim([-1 1]); zlim([-1 1])


myPolarPlot(thetas_obs,phis_obs,Edef_tot,2,clrmp,...
	'|E| = f(\Theta_{obs}-\Theta_{dip}) varying \Phi_{obs}')

myPolarPlot(thetas_obs,phis_obs,Edef_tot,1,clrmp,...
	'|E| = f(\Phi_{obs}) varying \Theta_{obs}')

myPolarPlot(thetas_obs,phis_obs,Edef_s,2,clrmp,...
	'|E_{s\Phi}| = f(\Theta_{obs})')

myPolarPlot(thetas_obs,phis_obs,Edef_s,1,clrmp,...
	'|E_{s\Phi}| = f(\Phi_{obs})')

nexttile
funprep

myPolarPlot(thetas_obs,phis_obs,Edef_p,2,clrmp,...
	'|E_{p\Theta}| = f(\Theta_{obs})')

myPolarPlot(thetas_obs,phis_obs,Edef_p,1,clrmp,...
	'|E_{p, u_\Theta}| = f(\Phi_{obs})')

sgtitle(sprintf('Simulation of E from dipole with \\Theta = %.0f° and \\Phi = %.0f°: %s',round(rad2deg(theta_dip)), round(rad2deg(phi_dip)),mode),'FontWeight','bold')
end


function myPolarPlot(theta, phi, E,n,clrmp, mytitle)
nexttile
for i=1:5:size(E,n)
	switch n
		case 2
			polarplot(theta,E(:,i));
		case 1
			polarplot(phi,E(i,:));
	end
	if i==1; hold on; end
end
switch n
	case 2
		thetalim([0 180]); title(mytitle);
	case 1
		title(mytitle)
end
ax=gca;
ll = length(ax.Children); cg = clrmp(ll);
ax.ColorOrder = cg;
end


function DOP = getDOP(x)
res = (max(x)-min(x))./(max(x)+min(x));
DOP = res;
end


function [x, y, z] = spherical2cart(r, theta, phi)
x = zeros(size(r,1),size(r,2));
y = zeros(size(r,1),size(r,2));
z = zeros(size(r,1),size(r,2));

for i=1:length(theta)
	for j=1:length(phi)
		% 		size(z)
		% 		r(i,j)
		% 		sin(theta(i))
		% 		sin(phi(j))
		x(i,j) = r(i,j)*sin(theta(i))*sin(phi(j));
		y(i,j) = r(i,j)*sin(theta(i))*cos(phi(j));
		z(i,j) = r(i,j)*cos(theta(i));
	end
end
end