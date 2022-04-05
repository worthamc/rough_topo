% noticed that rossby wave phase speeds from altimetry appear to change
% over time, and want to see if changes in hydrography could be
% responsible.  I have computed mean hydrography from ECCO2 over 3 year
% ranges.  Use that to compute modes at a particular point.
%
% Uses local average hydrography from OCCA or ECCO2 and ECCO2 topography
% spectrum.
%
% Should try again with HydroBase: http://www.whoi.edu/science/PO/hydrobase/index.html
%
% Based on the BPD_pde.m function.
%
% C. Wortham, 1 Sept 2020

% set this depending on where script is running
basedir = '/Users/worthamc/Documents/research/';
%basedir = '/data/worthamc/research/';
%addpath(genpath([basedir 'matlab/']),'-end')

% global stuff
global lon_qd lat_qd topo_x topo_y Dep_ll lon_hydro lat_hydro hydro_x hydro_y S_mean T_mean rC

%% settings and load reusables

% lat-lon for topography and coriolis parameter
lat_map = 30;
lon_map = 190;

% averaging area for topography.  This is the 'radius' of the box in deg.
topo_x = 2.5;
topo_y = 2.5;

% averaging area for hydrography.  This is the 'radius' of the box in deg.
hydro_x = 2.5;
hydro_y = 2.5;

% load 1/4 deg lat-lon grid from ECCO2
load([basedir 'scripts/spectrum/ecco2/cs510_runs/cube84/LatLon_qd.mat'])
lon_hydro = lon_qd;
lat_hydro = lat_qd;

% load model bathymetry from ECCO2
load([basedir 'data/ecco2/cs510_runs/cube84/grid/ll4/Dep_ll.mat']);
Dep_ll = circshift(Dep_ll,[180*4 0 0]); % shift longitude to run 0--360

% get ECCO2 vertical grid
load([basedir 'data/ecco2/cs510_runs/cube84/grid/cs510/dr.txt'])
rF=[0 -cumsum(dr)']'; % layer interace depths
rC=[ -cumsum(dr)']'+dr/2; % middle of layer depths
hydrosource = 'ecco2 3-year means';

tic
% get file names to loop over
datadirS = [basedir 'data/ecco2/cs510_runs/cube84/SALTanom/SALTanom_mean/'];
datadirT = [basedir 'data/ecco2/cs510_runs/cube84/THETA/THETA_mean/'];
FileNamesS = dir([datadirS 'SALTanom_mean_*.mat']);
FileNamesT = dir([datadirT 'THETA_mean_*.mat']);
% convert to cell
FileNamesS = {FileNamesS.name};
FileNamesT = {FileNamesT.name};
% setup some variables to be filled
omega_BPDtime = nan(1,length(FileNamesS));
omega_BTTtime = nan(1,length(FileNamesS));
k_time = nan(1,length(FileNamesS));
years = nan(2,length(FileNamesS));
N2cph_time = nan(length(rC),length(FileNamesS));
years = nan(2,length(FileNamesS));
for ii = 1:length(FileNamesS)
	load([datadirS FileNamesS{ii}])
	load([datadirT FileNamesT{ii}])
	
% if deeper than 1000 m ( rC(29)=-1007m )... and no land in the area... and
% not at the equator.
if ~isnan( T_mean(lon_hydro==lon_map , lat_hydro==lat_map , find(rC<-1000,1,'first') ))...
        && ~any(any( isnan(T_mean( abs(lon_hydro-lon_map)<=hydro_x , abs(lat_hydro-lat_map)<=hydro_y , 1 )) ))...
        && ~any(any( Dep_ll(abs(lon_qd-lon_map)<=topo_x , abs(lat_qd-lat_map)<=topo_y)==0 ))...
        && abs(lat_map) > 3

	% setting up a while loop to allow repeating the whole thing if
	% topography realization doesn't permit the desired mode. 
	while true
	
	disp(['ii = ' num2str(ii) ', lat = ' num2str(lat_map) ', lon = ' num2str(lon_map)])
	
	% average topography profile
% 	[h,y] = BPD_gettopo(lat_map,lon_map);
	h = H;
	
	% load hydrography and make average hydrography
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
    target_vals = 86400./omega_BTT(2); % converting the BTT frequency to equivalent eigenvalue alpha
    num_vals = 20; % number of eigenvalues to find
	[eigenvecs_mean,eigenvecs2D,omega_BPD,eigenvals,indBC] = BPD_sortmodes(S_temp',T_temp',rC_temp',lat_map,k,y,h,target_vals,num_vals,1);
	
	% if the output from BPD_sortmodes.m was empty, re-run with new topo
	if isempty(eigenvecs_mean)
		continue
	end
	
	break % breaking the while loop to continue to next iteration of for
	end
	
	% calculate the phase speed
	BPD_phasespeed = omega_BPD(indBC)/k; % this is in km/day
	disp(['phasespeed is ' num2str(BPD_phasespeed) ' km/day'])
% 	% and other output
	omega_BPDtime(ii) = omega_BPD(indBC);
	omega_BTTtime(ii) = omega_BTT(2);
	k_time(ii) = k;
	N2cph_time(1:length(rC_temp)-1,ii) = N2cphB;
	years(:,ii) = MeanYears;
	
end

toc
end


%% plot

figure(1);clf
plot(mean(years,1),-omega_BPDtime./k_time,'ko',mean(years,1),-omega_BTTtime./k_time,'k+')
xlabel('year')
ylabel('westward phase speed (km/day)')
legend('with topography','flat bottom')

figure(2);clf
plot(mean(years,1),omega_BPDtime./omega_BTTtime,'k.','markersize',15)
ylabel('c_{topo}/c_{flat}')
xlabel('year')

figure(3);clf
plot(sqrt(N2cph_time),rC)
xlabel('N (cph)')
ylabel('depth (m)')

figure(4);clf
plot(mean(years,1),nanmean(sqrt(N2cph_time),1))
ylabel('vertical average N (cph)')
xlabel('year')

figure(5);clf
plot(mean(years,1),nanmean(sqrt(N2cph_time(15:25,:)),1))
ylabel('thermocline average N (cph)')
xlabel('year')