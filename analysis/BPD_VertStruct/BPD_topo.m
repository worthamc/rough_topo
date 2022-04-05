% plots topography amplidude dependence of the BPD (small scale topography
% theory) 
%
% Uses local average hydrography from OCCA or ECCO2 and ECCO2 topography
% spectrum.
%
% Based on the BPD_pde.m function.
%
% C. Wortham, 19 May 2020

% set this depending on where script is running
basedir = '/Users/worthamc/Documents/research/';
%basedir = '/data/worthamc/research/';
%addpath(genpath([basedir 'matlab/']),'-end')
load([basedir 'matlab/CW/colormap/standard_cmap.mat'],'standard_cmap')
load([basedir 'matlab/CW/colormap/speed_cmap.mat'],'speed_cmap')
load([basedir 'matlab/CW/colormap/anomaly.mat'],'anomaly')

% global stuff
global lon_qd lat_qd topo_x topo_y Dep_ll lon_hydro lat_hydro hydro_x hydro_y S_mean T_mean rC

%% settings and load reusables

% lat-lon for topography and stratification
lat_map = 10;
lon_map = 190;

% averaging area for topography.  This is the 'radius' of the box in deg.
topo_x = 2.5;
topo_y = 2.5;

% averaging area for hydrography.  This is the 'radius' of the box in deg.
hydro_x = 2.5;
hydro_y = 2.5;

% load 1/4 deg lat-lon grid from ECCO2
load([basedir 'scripts/spectrum/ecco2/cs510_runs/cube84/LatLon_qd.mat'])

% load model bathymetry from ECCO2
load([basedir 'data/ecco2/cs510_runs/cube84/grid/ll4/Dep_ll.mat']);
Dep_ll = circshift(Dep_ll,[180*4 0 0]); % shift longitude to run 0--360

% get OCCA hydrography
% load OCCA grid (check locations!)
lat_hydro = nc_varget([basedir 'data/OCCA_1x1_v2/annual/DDtheta.0406annclim.nc'],'Latitude_t');
lat_hydro = lat_hydro-.5;
lon_hydro = nc_varget([basedir 'data/OCCA_1x1_v2/annual/DDtheta.0406annclim.nc'],'Longitude_t');
lon_hydro = lon_hydro-.5;
rC = -nc_varget([basedir 'data/OCCA_1x1_v2/annual/DDtheta.0406annclim.nc'],'Depth_c');
% load OCCA data
S_mean = nc_varget([basedir 'data/OCCA_1x1_v2/annual/DDsalt.0406annclim.nc'],'salt');
T_mean = nc_varget([basedir 'data/OCCA_1x1_v2/annual/DDtheta.0406annclim.nc'],'theta');
% rearrange dimensions to match ECCO2 style
S_mean = permute(S_mean,[3,2,1]);
T_mean = permute(T_mean,[3,2,1]);
% make OCCA stuff double
lat_hydro = double(lat_hydro);
lon_hydro = double(lon_hydro);
rC = double(rC);
S_mean = double(S_mean);
T_mean = double(T_mean);
hydrosource = 'occa';


%% new approach... 
%calculate phase speed for stratification and topography spectrum at a
%fixed point, but topography amplitude taken along a meridional line.  

tic
% loop for plotting dependence on topography amplitude.
load([basedir 'scripts/BPD_VertStruct/BPD_phasespeed_global.mat'],'h_rms')
h_rms = nanmean(h_rms,2);
dLatLon = 2;
lat_topo = -80:dLatLon:78;
omega_BPDtopo = nan(1,length(h_rms));
omega_BPDtopo2 = nan(1,length(h_rms));
omega_BTTtopo = nan(1,length(h_rms));
k_topo = nan(1,length(h_rms));
for ii = 1:length(h_rms)
	if ~isnan(h_rms(ii))
    amp = h_rms(ii);

% if deeper than 1000 m ( rC(29)=-1007m )... and no land in the area... and
% not at the equator.
if ~isnan( T_mean(lon_hydro==lon_map , lat_hydro==lat_map , find(rC<-1000,1,'first') ))...
        && ~any(any( isnan(T_mean( abs(lon_hydro-lon_map)<=hydro_x , abs(lat_hydro-lat_map)<=hydro_y , 1 )) ))...
        && ~any(any( Dep_ll(abs(lon_qd-lon_map)<=topo_x , abs(lat_qd-lat_map)<=topo_y)==0 ))...
        && abs(lat_map) >= 3

	% setting up a while loop to allow repeating the whole thing if
	% topography realization doesn't permit the desired mode. 
	while true
	
	disp(['ii = ' num2str(ii) ', lat = ' num2str(lat_map) ', lon = ' num2str(lon_map)])
	
	% average topography profile
% 	[h,y] = BPD_gettopo(lat_map,lon_map);
	h = amp*(H-mean(H))/sqrt(var(H))+ mean(H);
	
	% average hydrography
	[S_temp,T_temp,rC_temp] = BPD_gethydro(lat_map,lon_map);
	
	% plot BTT modes, and use that for wavelength
    path(path,[basedir 'scripts/VertStruct']);
    [LdB,eigenvecsB,rhoB,N2cphB] = BTTmodes(S_temp',T_temp',rC_temp',lat_map,1,5); % second to last argument zero or one for plotting
    k = 1/LdB(2); % (cpk) 
    k = -k/50 % something like -1/10x the deformation wavenumber
   	% the BTT frequencies
	BETA = 2*7.292e-5*cosd(lat_map)/6370e3;
    omega_BTT = -BETA*(k*2*pi/1000)./((k*2*pi/1000)^2 + 1./(LdB*1000).^2)*86400/2/pi; % cycles/day

	% solve the problem, looking for solutions with eigenvalue near the BTT
    % first baroclinic mode frequency
    target_vals = omega_BTT(2); % will search for eigenvectors with eigenvalues close to the first baroclinic mode frequency
    num_vals = 20; % number of eigenvalues to find
	[eigenvecs_mean,eigenvecs2D,w_r,w_i,eigenvals,indBC] = BPD_sortmodes(S_temp',T_temp',rC_temp',lat_map,k,y,h,target_vals,num_vals);
	
	% if the output from BPD_sortmodes.m was empty, re-run with new topo
	if isempty(eigenvecs_mean)
		continue
	end
	
	break % breaking the while loop to continue to next iteration of for
	end
	
	% calculate the phase speed
	BPD_phasespeed = w_r(indBC)/k; % this is in km/day
	disp(['phasespeed is ' num2str(BPD_phasespeed) ' km/day'])
	% and other output
	omega_BPDtopo(ii) = w_r(indBC);
	omega_BPDtopo2(ii) = w_r(7);
	omega_BTTtopo(ii) = omega_BTT(2);
	k_topo(ii) = k;
	
end
	
end

toc
end


% plot
figure(1);clf
plot(lat_topo,omega_BPDtopo/omega_BTT(2),'k.','markersize',15)
xlabel('latitude of topography amplitude')
ylabel('c_{topo}/c_{flat}')
xlim([-60 60])
ylim([1 2])

figure(2);clf
plot(lat_topo,-omega_BPDtopo./k_topo,'ko',lat_topo,-omega_BTTtopo./k_topo,'k+')
xlabel('year')
ylabel('westward phase speed (km/day)')
legend('with topography','flat bottom')


return

%% old approach... 
% just set amplitude of topography arbitrarily.  use fixed topography
% spectrum and stratification.  

tic
% loop for plotting dependence on topography amplitude.
AMPmax = 1000;
Namp = 30;
omega_BPDtopo = nan(1,Namp);
omega_BPDtopo2 = nan(1,Namp);
omega_BTTtopo = nan(1,Namp);
k_topo = nan(1,Namp);
AMP = linspace(10,AMPmax,Namp); % multiplicative factor for the amplitude
for ii = 1:Namp
    amp = AMP(ii);

% if deeper than 1000 m ( rC(29)=-1007m )... and no land in the area... and
% not at the equator.
if ~isnan( T_mean(lon_hydro==lon_map , lat_hydro==lat_map , find(rC<-1000,1,'first') ))...
        && ~any(any( isnan(T_mean( abs(lon_hydro-lon_map)<=hydro_x , abs(lat_hydro-lat_map)<=hydro_y , 1 )) ))...
        && ~any(any( Dep_ll(abs(lon_qd-lon_map)<=topo_x , abs(lat_qd-lat_map)<=topo_y)==0 ))...
        && abs(lat_map) >= 3

	% setting up a while loop to allow repeating the whole thing if
	% topography realization doesn't permit the desired mode. 
	while true
	
	disp(['ii = ' num2str(ii) ', lat = ' num2str(lat_map) ', lon = ' num2str(lon_map)])
	
	% average topography profile
% 	[h,y] = BPD_gettopo(lat_map,lon_map);
% 	H0 = mean(h);
% 	H = h - H0;
% 	H = H/max(H); % save normalized H
% 	h = H0 + amp*H;

	% make some simple topography
	h = -5000 + amp*cos(2*pi*4*y/max(y));
	
	% average hydrography
	[S_temp,T_temp,rC_temp] = BPD_gethydro(lat_map,lon_map);
	
	% plot BTT modes, and use that for wavelength
    path(path,[basedir 'scripts/VertStruct']);
    [LdB,eigenvecsB,rhoB,N2cphB] = BTTmodes(S_temp',T_temp',rC_temp',lat_map,1,5); % second to last argument zero or one for plotting
    k = 1/LdB(2); % (cpk) 
    k = -k/50 % something like -1/10x the deformation wavenumber
   	% the BTT frequencies
	BETA = 2*7.292e-5*cosd(lat_map)/6370e3;
    omega_BTT = -BETA*(k*2*pi/1000)./((k*2*pi/1000)^2 + 1./(LdB*1000).^2)*86400/2/pi; % cycles/day

	% solve the problem, looking for solutions with eigenvalue near the BTT
    % first baroclinic mode frequency
    target_vals = omega_BTT(2); % will search for eigenvectors with eigenvalues close to the first baroclinic mode frequency
    num_vals = 100; % number of eigenvalues to find
	[eigenvecs_mean,eigenvecs2D,w_r,w_i,eigenvals,indBC] = BPD_sortmodes(S_temp',T_temp',rC_temp',lat_map,k,y,h,target_vals,num_vals,1);
	
	% if the output from BPD_sortmodes.m was empty, re-run with new topo
	if isempty(eigenvecs_mean)
		continue
	end
	
	break % breaking the while loop to continue to next iteration of for
	end
	
	% identify the modes with y-dependence (change of sign)
	ymode = ~(abs(sum(squeeze(sign(real(eigenvecs2D(1,:,:)))),1))==size(eigenvecs2D(1,:,:),2));
	% identify modes with n interior extrema 
	X = sign(diff(real(eigenvecs_mean),1,1)); % this is a rough derivative
	Y = X(1:end-1,:) + X(2:end,:); % Y=0 when the derivative changes sign... at extrema
	NumMax = sum(Y==0,1); % counting the number of interior extrema
	% now save the indices of barotropic and baroclinic modes with least meridional structure.
	VertMode = [1]; % start with BT mode
	NumMax(1) = -1; % exclude BT mode from further consideration
	for jj=0:max(NumMax)
		ind_temp = find(NumMax==jj & ~ymode);
%  		% pick mode with largest range
% 		[~,index] = max(range(eigenvecs_mean(:,NumMax==jj),1));
% 		% pick mode with largest surface amplitude
% 		[~,index] = max(eigenvecs_mean(1,NumMax==jj));
		% pick mode with least variance at bottom
		[~,index] = min(var(squeeze(eigenvecs2D(end,:,ind_temp)),1)); % variance
% 		[~,index] = min(sqrt(mean(squeeze(eigenvecs2D(end,:,ind_temp).^2)))); % rms
		VertMode = [VertMode ind_temp(index)];
	end
  	indBC = VertMode(2);
	
	% calculate the phase speed
	BPD_phasespeed = w_r(indBC)/k; % this is in km/day
	disp(['phasespeed is ' num2str(BPD_phasespeed) ' km/day'])
	% and other output
	omega_BPDtopo(ii) = w_r(indBC);
	omega_BPDtopo2(ii) = w_r(7);
	omega_BTTtopo(ii) = omega_BTT(2);
	k_topo(ii) = k;
	
	ind_save(ii) = indBC;
	wr_save(ii,:) = w_r;
		
end

toc
end



% plot

% % use this plotting code at first to plot histogram of h_rms.
% load([basedir 'scripts/BPD_VertStruct/BPD_phasespeed_global.mat'],'h_rms')
% figure(1);clf
% hold on
% % multiplying h_rms by sqrt(2) for direct comparison with amplitude of
% % sinusoidal topography (alternatively could divide topography amplitude by
% % sqrt(2).
% [n,xout] = hist(h_rms(:)*sqrt(2),100); 
% b1 = bar(xout,n); % plot histogram for h_rms
% ax1 = gca;
% set(ax1,'xaxislocation','top','yaxislocation','right','xticklabel',[]) % set axes position
% ylabel('\# of grid points')
% box off
% ax2 = axes('position',get(ax1,'position'),'yaxislocation','left','color','none'); % make second axes
% line(AMP,omega_BPDtopo/omega_BTT(2),'color','k','linewidth',2)
% text(.9*AMP(end),omega_BPDtopo(end)/omega_BTT(2),[num2str(lat_map) '^\circ'])
% set(ax1,'xlim',get(ax2,'xlim')) % make x axes match
% set(ax1,'position',get(ax2,'position')) % and really make sure they match!
% set(b1,'edgecolor',[.5 .5 .5],'facecolor',[.5 .5 .5]) % set color of histogram
% xlabel('topography amplitude (m)')
% ylabel('c_{topo}/c_{flat}')

% use this plotting code to add more lines.
figure(1);
hold on
plot(AMP,omega_BPDtopo/omega_BTT(2),'k','linewidth',2)
xlabel('topography amplitude (m)')
ylabel('c_{topo}/c_{flat}')
text(.9*AMP(end),omega_BPDtopo(end)/omega_BTT(2),[num2str(lat_map) '^\circ'])
% text(AMP(end),omega_BPDtopo(end)/omega_BTT(2),[num2str(2) '^\circ'])


figure(2);
% clf
hold on
plot(-AMP/mean(h),omega_BPDtopo/omega_BTT(2),'k')
xlabel('fractional amplitude')
ylabel('c_{topo}/c_{flat}')

figure(3);
% clf
hold on
plot(AMP,-omega_BPDtopo./k_topo,'ko',AMP,-omega_BTTtopo./k_topo,'k+')
xlabel('topography amplitude (m)')
ylabel('westward phase speed (km/day)')
legend('with topography','flat bottom')

%% different approach, 
% using maps of phase speed with topography amplitude halved or doubled
% These maps are saved as BPD_phasespeed_05xTOPO.mat,
% BPD_phasespeed_global.mat, and BPD_phasespeed_2xTOPO.mat.  

% load the save calculations
load([basedir 'scripts/BPD_VertStruct/BPD_phasespeed_05xTOPO.mat'])
BPD_phasespeed_05x = BPD_phasespeed;
BPD_omega_05x = BPD_omega;
load([basedir 'scripts/BPD_VertStruct/BPD_phasespeed_global.mat'])
BPD_phasespeed_1x = BPD_phasespeed;
BPD_omega_1x = BPD_omega;
load([basedir 'scripts/BPD_VertStruct/BPD_phasespeed_2xTOPO.mat'])
BPD_phasespeed_2x = BPD_phasespeed;
BPD_omega_2x = BPD_omega;

load([basedir 'matlab/CW/colormap/anomaly.mat'],'anomaly') % load my anomaly colormap

figure(1);clf
m_proj('miller cylindrical','lon',[0 360],'lat',[-80 80]); % or robinson
m_pcolor(lon_map-dLatLon/2,lat_map-dLatLon/2,BPD_phasespeed_05x./BPD_phasespeed_1x);shading flat;
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in');
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('halved topography')
colorbar
set(gcf,'colormap',anomaly)
caxis([.7 1.3])

figure(2);clf
m_proj('miller cylindrical','lon',[0 360],'lat',[-80 80]); % or robinson
% m_pcolor(lon_map-dLatLon/2,lat_map-dLatLon/2,BPD_phasespeed_2x./BPD_phasespeed_1x);shading flat;
m_contourf(lon_map-dLatLon/2,lat_map-dLatLon/2,BPD_phasespeed_2x./BPD_phasespeed_1x,20,'linecolor','none')
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in');
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('doubled topography')
colorbar
set(gcf,'colormap',anomaly)
caxis([.7 1.3])


BTT_phasespeed = BTT_omega./wavenumber;

% ratio
ratio = BPD_phasespeed_2x./BTT_phasespeed;
ratio(ratio>4.5)=4.5; 
figure(3);clf
m_proj('miller cylindrical','lon',[0 360],'lat',[-80 80]); % or robinson
m_contourf(lon_map-dLatLon/2,lat_map-dLatLon/2,ratio,20,'linecolor','none');
caxis([1 4.5])
plot_axes = gca;
cbar_axes = colorbar;
% c_map = colormap(jet(20));
c_map = colormap(standard_cmap);
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in');
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('phase speed ratio')
cbar_triangles(cbar_axes,plot_axes,c_map,1,0,0,0)

ratio = BPD_phasespeed_05x./BTT_phasespeed;
ratio(ratio>4.5)=4.5; 
figure(4);clf
m_proj('miller cylindrical','lon',[0 360],'lat',[-80 80]); % or robinson
m_contourf(lon_map-dLatLon/2,lat_map-dLatLon/2,ratio,20,'linecolor','none');
caxis([1 4.5])
plot_axes = gca;
cbar_axes = colorbar;
% c_map = colormap(jet(20));
c_map = colormap(standard_cmap);
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in');
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('phase speed ratio')
cbar_triangles(cbar_axes,plot_axes,c_map,1,0,0,0)

