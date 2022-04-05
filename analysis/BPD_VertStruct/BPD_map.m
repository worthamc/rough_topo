% plots a map of the BPD (small scale topography theory) long-wave phase
% speed.
% Compare to zonal phase speed propagation in the model or to Fu (2009)
% Fig. 11
%
% Uses local average hydrography from OCCA or ECCO2 and ECCO2 topography
% spectrum.
%
% the older version (BPD_pde_map.m) tried to do this with ECCO2 topography
% and hydrography, but there were too many problems with he salinity
% fields.
%
% Would like to add option to use Smith and Sandwell topography, as well.
%
% Based on the BPD_pde.m function.
%
% C. Wortham, 21 Aug. 2020

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

% lat-lon grid for map... lower left corners
dLatLon = 2; 
lat_map = -80:dLatLon:78;
lon_map = 0:dLatLon:358;
[lon_map,lat_map] = meshgrid(lon_map,lat_map);

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

% initialize
BPD_phasespeed = nan(size(lat_map)); % save phase speed
BPD_mode = nan([size(lat_map),50]); % save vertical structures
BTT_mode1 = nan([size(lat_map),50]);
BPD_omega = nan(size(lat_map)); % save frequencies
BTT_omega = nan(size(lat_map));
wavenumber = nan(size(lat_map)); % save wavenumber
h_rms = nan(size(lat_map)); % save rms of topography, as proxy for bottom roughness

tic
% for ii=1:length(lat_map(:)) % using linear indexing
% for ii=find(lon_map>=270.1 & lon_map<=360)'
% for ii=find(lon_map>=120.1 & lon_map<=300)'
% for ii=find(lon_map==196)'
% for ii=find(lat_map==30 & lon_map==210) 
for ii=find(lat_map==16 & lon_map==210) 

	
% if deeper than 1000 m ( rC(29)=-1007m )... and no land in the area... and
% not at the equator.
if ~isnan( T_mean(lon_hydro==lon_map(ii) , lat_hydro==lat_map(ii) , find(rC<-1000,1,'first') ))...
        && ~any(any( isnan(T_mean( abs(lon_hydro-lon_map(ii))<=hydro_x , abs(lat_hydro-lat_map(ii))<=hydro_y , 1 )) ))...
        && ~any(any( Dep_ll(abs(lon_qd-lon_map(ii))<=topo_x , abs(lat_qd-lat_map(ii))<=topo_y)==0 ))...
        && abs(lat_map(ii)) >= 3

	% setting up a while loop to allow repeating the whole thing if
	% topography realization doesn't permit the desired mode. 
	while true
	
	disp(['ii = ' num2str(ii) ', lat = ' num2str(lat_map(ii)) ', lon = ' num2str(lon_map(ii))])
	
	% average topography profile
	[h,y] = BPD_gettopo(lat_map(ii),lon_map(ii));
	% using saved topography
% 	h=H;
	% make some simple topography
% 	h = mean(h) + sqrt(var(h))*cos(2*pi*4*y/max(y));
% 	h = mean(h) - sqrt(var(h))*linspace(-1,1,length(h)); % negative slope... won't find BC solution faster than BTT BC1 mode!
% 	h = mean(h) + sqrt(var(h))*linspace(-1,1,length(h)); % positive slope


	% average hydrography
	[S_temp,T_temp,rC_temp] = BPD_gethydro(lat_map(ii),lon_map(ii));
	% using saved hydrography
% 	S_temp = S';
% 	T_temp = T';
% 	rC_temp = rC;
	% trying to creating barotropic hydrographic profile... but failing
% 	S_temp = mean(S_temp)*ones(size(S_temp)); % just use average salinity
% 	pres = sw_pres(-rC_temp,lat_map(ii));
% 	T_temp = sw_temp(S_temp,mean(T_temp)*ones(size(S_temp)),1.1*pres,0); clear pres;% and an in situ temp profile that will yield constant potential temp
% 	T_temp = mean(T_temp)*ones(size(S_temp));


	% plot BTT modes, and use that for wavelength
    path(path,[basedir 'scripts/VertStruct']);
    [LdB,eigenvecsB,rhoB,N2cphB] = BTTmodes(S_temp',T_temp',rC_temp',lat_map(ii),1,5); % second to last argument zero or one for plotting
    k = 1/LdB(2); % (cpk) 
    k = -k/50 % something like -1/10x the deformation wavenumber
    % the BTT frequencies
	BETA = 2*7.292e-5*cosd(lat_map(ii))/6370e3;
    omega_BTT = -BETA*(k*2*pi/1000)./((k*2*pi/1000)^2 + 1./(LdB*1000).^2)*86400/2/pi; % cycles/day

	% solve the problem, looking for solutions with eigenvalue near the BTT
    % first baroclinic mode frequency
    target_vals = omega_BTT(2); % will search for eigenvectors with eigenvalues close to the first baroclinic mode frequency
    num_vals = 100; % number of eigenvalues to find
	[eigenvecs_mean,eigenvecs2D,w_r,w_i,eigenvals,indBC] = BPD_sortmodes(S_temp',T_temp',rC_temp',lat_map(ii),k,y,h,target_vals,num_vals,1);
	
	% if the output from BPD_sortmodes.m was empty, re-run with new topo
	if isempty(eigenvecs_mean)
		continue
	end
	
	break % breaking the while loop to continue to next iteration of for
	end
	
% 	% identify modes with n interior extrema 
% 	X = sign(diff(real(eigenvecs_mean),1,1)); % this is a rough derivative
% 	Y = X(1:end-1,:) + X(2:end,:); % Y=0 when the derivative changes sign... at extrema
% 	NumMax = sum(Y==0,1); % counting the number of interior extrema
% 	% now save the indices of barotropic and baroclinic modes with least meridional structure.
% 	VertMode = [1]; % start with BT mode
% 	for ii=0:max(NumMax)
% 		ind_temp = find(NumMax==ii);
% 		% pick mode with largest amplitude
% 		[~,index] = max(range(eigenvecs_mean(:,NumMax==ii),1));
% 		VertMode = [VertMode ind_temp(index)];
% 	end

% 	% identify modes with n zero crossings
% 	X = sign(real(eigenvecs_mean)); % sign of function
% 	Y = X(1:end-1,:) + X(2:end,:); % count changes in sign... equal to number of zero crossings
% 	NumZero = sum(Y==0,1); % counting the number of zero crossings
% 	% now save the indices of barotropic and baroclinic modes with least meridional structure.
% 	VertMode = 1; % start with BT and gravest baroclinic modes
% 	for ii=1:max(NumZero)
% 		ind_temp = find(NumZero==ii);
% 		% pick mode with largest amplitude?
% 		[~,index] = max(range(eigenvecs_mean(:,NumZero==ii),1));
% 		VertMode = [VertMode ind_temp(index)];
% 	end

	% calculate the phase speed
	BPD_phasespeed(ii) = w_r(indBC)/k; % this is in km/day
	disp(['phasespeed is ' num2str(BPD_phasespeed(ii)) ' km/day'])

	% save some other output
	[a,b] = ind2sub(size(lat_map),ii);
	BPD_mode(a,b,1:length(rC_temp)) = eigenvecs_mean(:,indBC);
	BTT_mode1(a,b,1:length(rC_temp)) = eigenvecsB(:,2);
	BPD_omega(ii) = w_r(indBC);
	BTT_omega(ii) = omega_BTT(2);
	wavenumber(ii) = k;
	h_rms(ii) = sqrt(mean((h-mean(h)).^2));
	
end

toc
end

save('~/Desktop/BPD_phasespeed.mat','BPD_phasespeed','BPD_mode','BTT_mode1','BPD_omega','BTT_omega','wavenumber','h_rms','lat_map','lon_map','dLatLon','hydrosource','hydro_x','hydro_y','topo_x','topo_y')

%% draw map

figure(1);clf
m_proj('miller cylindrical','lon',[0 360],'lat',[-80 80]); % or robinson
% m_pcolor(lon_map,lat_map,log10(-BPD_phasespeed));shading flat;
m_pcolor(lon_map-dLatLon/2,lat_map-dLatLon/2,BPD_phasespeed);shading flat;
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in');
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('BPD phase speed (km/day)')
colorbar
colormap(jet(20))
caxis([-10 0])


return
%% some plan-view plots of the solution

fno = 30;
% y is in meters
% k is in cpk
% x is in meters
% z is in meters

% the solution... select 1 for BT mode, or indBC for baroclinic mode
solnYZ = eigenvecs2D(:,:,1)';

% the wave part
x = 1000*linspace(0,-1/k,50);
waveX = -sin(2*pi*k/1000*x)';

% make them 3D matrices and full solution
waveX = repmat(waveX,[1,size(solnYZ,1),size(solnYZ,2)]);
solnYZ = repmat(solnYZ,[1,1,size(waveX,1)]);
solnYZ = permute(solnYZ,[3,1,2]);
solnXYZ = solnYZ.*waveX;

% the topography
topoXY = repmat(h,[size(waveX,1),1]);

% the velocity
dx = x(2)-x(1);
dy = y(2)-y(1);
F = sw_f(lat_map(ii));
RHO0 = 1000;
v = 1/RHO0/F * diff(solnXYZ,1,1)/dx;
u = -1/RHO0/F * diff(solnXYZ,1,2)/dy;
% V = 1/RHO0/F * diff(waveX,1,1)/dx;
% U = -1/RHO0/F * diff(waveX,1,2)/dy;
V = mean(v,2);
V = repmat(V,[1 size(v,2) 1]);
U = mean(u,2);
U = repmat(U,[1 size(u,2) 1]);

fno = fno+1;
figure(fno);clf;
pcolor(x,y,waveX(:,:,1)')
title('the large scale wave')
colorbar

fno = fno+1;
figure(fno);clf;
surf(x,y,topoXY');shading flat
title('the topography')
colorbar
colormap(bone)

iplot = 1; % select index of depth to plot

fno = fno+1;
figure(fno);clf;
pcolor(x,y,squeeze(solnYZ(:,:,iplot))')
title('the 2D solution at z= ??')
colorbar
caxis([-2.5 2.5])

fno = fno+1;
figure(fno);clf;
pcolor(x,y,squeeze(solnXYZ(:,:,iplot)'));shading flat
title('full solution at z= ??')
colorbar
caxis([-2.5 2.5])
hold on
quiver(x(1:end-1),y(1:end-1),squeeze(u(1:end-1,:,iplot))',squeeze(v(:,1:end-1,iplot))',2)
hold off


fno = fno+1;
figure(fno);clf;
subplot(1,4,1)
plot(h,y)
subplot(1,4,[2 4])
pcolor(x,y,squeeze(solnXYZ(:,:,iplot)'));shading flat
title('"induced flow" at z= ??')
colorbar
caxis([-2.5 2.5])
hold on
quiver(x(1:end-1),y(1:end-1),squeeze(u(1:end-1,:,iplot)-U(1:end-1,:,iplot))',squeeze(v(:,1:end-1,iplot)-V(:,1:end-1,iplot))',2)
hold off
