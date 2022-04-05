% plots stratification dependence of the BPD (small scale topography
% theory).  This computes the phase speed using fixed coriolis parameter
% and topography, but taking stratification from different latitudes.
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

% global stuff
global lon_qd lat_qd topo_x topo_y Dep_ll lon_hydro lat_hydro hydro_x hydro_y S_mean T_mean rC

%% settings and load reusables

% lat-lon for topography and coriolis parameter
% lat_map = 30;
% lon_map = 190;
lat_map = 26;
lon_map = 210;

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

% % use zonal average hydrography
% S_mean = repmat(nanmean(S_mean,1),[length(lon_hydro),1,1]);
% T_mean = repmat(nanmean(T_mean,1),[length(lon_hydro),1,1]);

tic
% loop for plotting dependence on stratification.
dLatLon = 2;
LAT_FAKE = -50:dLatLon:50;
% LAT_FAKE = 42;
omega_BPDhydro = nan(1,length(LAT_FAKE));
omega_BTThydro = nan(1,length(LAT_FAKE));
k_hydro = nan(1,length(LAT_FAKE));
for ii = 1:length(LAT_FAKE)
    lat_fake = LAT_FAKE(ii);

% if deeper than 1000 m ( rC(29)=-1007m )... and no land in the area... and
% not at the equator.
if ~isnan( T_mean(lon_hydro==lon_map , lat_hydro==lat_fake , find(rC<-1000,1,'first') ))...
        && ~any(any( isnan(T_mean( abs(lon_hydro-lon_map)<=hydro_x , abs(lat_hydro-lat_fake)<=hydro_y , 1 )) ))...
        && ~any(any( Dep_ll(abs(lon_qd-lon_map)<=topo_x , abs(lat_qd-lat_map)<=topo_y)==0 ))...
        && abs(lat_fake) > 3

	% setting up a while loop to allow repeating the whole thing if
	% topography realization doesn't permit the desired mode. 
	while true
	
	disp(['ii = ' num2str(ii) ', lat = ' num2str(lat_fake) ', lon = ' num2str(lon_map)])
	
	% average topography profile
% 	[h,y] = BPD_gettopo(lat_map,lon_map);
	h = H;
	
	% average hydrography... from the "fake" latitude
	[S_temp,T_temp,rC_temp] = BPD_gethydro(lat_fake,lon_map);
	
	% make the depth of the hydrography profiles 4000m
	if rC_temp(end)<-4000
		S_temp = S_temp(rC_temp>-4000);
		T_temp = T_temp(rC_temp>-4000);
		rC_temp = rC_temp(rC_temp>-4000);
	else
		S_temp = interp1(rC_temp,S_temp,rC(rC>-4000),'linear','extrap');
		T_temp = interp1(rC_temp,T_temp,rC(rC>-4000),'linear','extrap');
		rC_temp = rC(rC>-4000);
	end
	
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
	[eigenvecs_mean,eigenvecs2D,w_r,w_i,eigenvals,indBC] = BPD_sortmodes(S_temp',T_temp',rC_temp',lat_map,k,y,h,target_vals,num_vals,1);
	
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
	omega_BPDhydro(ii) = w_r(indBC);
	omega_BTThydro(ii) = omega_BTT(2);
	k_hydro(ii) = k;
	
end

toc
end


%% plot
figure(1);clf
plot(LAT_FAKE,omega_BPDhydro./omega_BTThydro,'k.','markersize',15)
xlabel('latitude of hydrographic profile')
ylabel('c_{topo}/c_{flat}')
xlim([-60 60])
ylim([1 2])

figure(2);clf
plot(LAT_FAKE,-omega_BPDhydro./k_hydro,'ko',LAT_FAKE,-omega_BTThydro./k_hydro,'k+')
xlabel('latitude of hydrographic profile')
ylabel('westward phase speed (km/day)')
legend('with topography','flat bottom')
xlim([-60 60])



