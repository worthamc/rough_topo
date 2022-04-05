% Similar to ECCO_4dspec.m but using the BPD vertical modes.  
% This takes 1/4deg gridded ECCO velocity profiles, computes the BPD
% vertical modes at each point, and decomposes the velocity profile into
% those modes.  Then compute come statistics.
%
% Right now, just doing the modal decomposition with first two modes.
% Getting pretty rediculous modal coefficients when I use 20 modes.  Really
% want to look at the vertical structure of the variance.  Do the BPD modes
% do a better job of capturing the verical structure of variance than the
% BTT modes?  Also, are the BPD modes any less coherent than the BTT first
% and second modes were?  These look like important tests of the theory.
%
% To get a better decomposition, I think I need to have more vertical
% modes.  That means more sophisticated sorting to find higher vertical
% modes without strong vertical structure.  Should at least be able to find
% n=2.
%
% C. Wortham, 19 May 2020

% set this depending on where script is running
% basedir = '/data/worthamc/research/';
basedir = '/Users/worthamc/Documents/research/';

% addpath(genpath([basedir 'matlab/']),'-end')

tic

lonvec = [180 220]; % east to west
latvec = [30 35]; % south to north
basin = 'Npac'; 

% import the data from the .data files, or load imported data
% [u,v,S,T,Dep_ll,lat,lon,time] = ECCO_uv_bin2mat(latvec,lonvec,4,4,1);
% SaveDate = date;
% save([basedir 'data/ecco2/cs510_runs/cube84/saved/data_uv_lat_' num2str(latvec(1)) '_' num2str(latvec(2)) '_lon_' num2str(lonvec(1)) '_' num2str(lonvec(2))],'u','v','S','T','Dep_ll','lat','lon','time','SaveDate');
filename = [basedir 'data/ecco2/cs510_runs/cube84/saved/data_uv_lat_' num2str(latvec(1)) '_' num2str(latvec(2)) '_lon_' num2str(lonvec(1)) '_' num2str(lonvec(2))];
load(filename);

toc

fno = 0;

% velocity has NaN for missing values

% % decrease resolution (no filtering)
% ny=4;
% nx=4;
% nt=1;
% u = u(1:ny:end,1:nx:end,:,1:nt:end);
% v = v(1:ny:end,1:nx:end,:,1:nt:end);
% S = S(1:ny:end,1:nx:end,:);
% T = T(1:ny:end,1:nx:end,:);
% Dep_ll = Dep_ll(1:ny:end,1:nx:end,:);
% time = time(1:nt:end);
% lon = lon(1:nx:end);
% lat = lat(1:ny:end);

% remove velocity time mean at each point
TimeMeanU = nanmean(u,4);
u = u - repmat(TimeMeanU, [1 1 1 size(u,4)]);
TimeMeanV = nanmean(v,4);
v = v - repmat(TimeMeanV, [1 1 1 size(v,4)]);

% some basic constants
Nreps = length(time);
Nlon = length(lon);
Nlat = length(lat);
dt = time(2) - time(1); % days
dx = lon(2)-lon(1); % deg
dy = lat(2) - lat(1); % deg
[distX,phase] = sw_dist([lat(floor(end/2)) lat(floor(end/2))],[0 1],'km');

% more needs... layer depths, bathymetry, mean S and T fields
load([basedir 'data/ecco2/cs510_runs/cube84/grid/cs510/dr.txt'])
rF=[0 -cumsum(dr)']'; % layer interace depths
rC=[ -cumsum(dr)']'+dr/2; % middle of layer depths
% make salt anomaly --> salt and potential temp --> in situ temp
S = S + 35;
pres = sw_pres(-rC,mean(lat)); pres = permute(pres,[3 2 1]); pres = repmat(pres,[length(lat) length(lon) 1]);
T = sw_temp(S,T,pres,zeros(size(pres))); clear pres;


% immitate mooring, with only 6 depths... need to uncomment this and
% a few other lines with "instrument" or "inst" for this to work.
% instruments = [11 22 25 34 42 45];


%% now get the stuff for BPD modes.  Copied from BPD_map.  Hope there are no conflicts.

%addpath(genpath([basedir 'matlab/']),'-end')
load([basedir 'matlab/CW/colormap/standard_cmap.mat'],'standard_cmap')
load([basedir 'matlab/CW/colormap/speed_cmap.mat'],'speed_cmap')
load([basedir 'matlab/CW/colormap/anomaly.mat'],'anomaly')

% global stuff
global lon_qd lat_qd topo_x topo_y Dep_ll lon_hydro lat_hydro hydro_x hydro_y S_mean T_mean rC

% settings and load reusables

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

%% compute vertical modes

path(path,[basedir 'scripts/VertStruct'])
ceil = 11; % only use data from below 100 meters for modal composition.
iplot = 1;
nplot = 5;
nmodes = 2;
alpha_u = nan(Nlat,Nlon,nmodes,Nreps); % modal composition.  Only computing first 10 modes.
alpha_v = nan(Nlat,Nlon,nmodes,Nreps); % modal composition.  Only computing first 10 modes.
for yy=1%:Nlat
	for xx=1:50%Nlon
		if ~isnan(u(yy,xx,nmodes,1)) && (Dep_ll(yy,xx) < -2000) % only look where you have deep velocities
			%(Dep_ll(yy,xx) < -2000); % only look where depth is more than 1000m 
			ind = find( ~isnan(squeeze(u(yy,xx,:,1))) ,1,'last' ); % index of last real velocity
			
			% pick out particular instruments
%			ind2 = find(instruments<ind,1,'last');
%			inst = instruments(1:ind2);
			

% compute BPD modes
	% average topography profile
	[h,y] = BPD_gettopo(lat(yy),lon(xx));
	% using saved topography
% 	h=H;
	% make some simple topography
% 	h = mean(h) + sqrt(var(h))*cos(2*pi*4*y/max(y));
% 	h = mean(h) - sqrt(var(h))*linspace(-1,1,length(h)); % negative slope... won't find BC solution faster than BTT BC1 mode!
% 	h = mean(h) + sqrt(var(h))*linspace(-1,1,length(h)); % positive slope


	% average hydrography
	[S_temp,T_temp,rC_temp] = BPD_gethydro(lat(yy),lon(xx));
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
    [LdB,eigenvecsB,rhoB,N2cphB] = BTTmodes(S_temp',T_temp',rC_temp',lat(yy),1,5,1); % second to last argument zero or one for plotting
    k = 1/LdB(2); % (cpk) 
    k = -k/50 % something like -1/10x the deformation wavenumber
    % the BTT frequencies
	BETA = 2*7.292e-5*cosd(lat_map(yy))/6370e3;
    omega_BTT = -BETA*(k*2*pi/1000)./((k*2*pi/1000)^2 + 1./(LdB*1000).^2)*86400/2/pi; % cycles/day

	% solve the problem, looking for solutions with eigenvalue near the BTT
    % first baroclinic mode frequency
    target_vals = omega_BTT(2); % will search for eigenvectors with eigenvalues close to the first baroclinic mode frequency
    num_vals = 20; % number of eigenvalues to find
	[eigenvecs_mean,eigenvecs2D,w_r,w_i,eigenvals,indBC] = BPD_sortmodes(S_temp',T_temp',rC_temp',lat(yy),k,y,h,target_vals,num_vals);
	
	% make the gravest baroclinic mode second
	eigenvecs = [ eigenvecs_mean(:,1) circshift(eigenvecs_mean(:,2:end),[0,-indBC+2]) ];

			% least squares fit	
			E_temp = eigenvecs(ceil:ind,1:nmodes);
			u_temp = squeeze(u(yy,xx,ceil:ind,:));
%			v_temp = squeeze(u(yy,xx,inst,:));			
			v_temp = squeeze(v(yy,xx,ceil:ind,:));			
%			v_temp = squeeze(v(yy,xx,inst,:));	
			alpha_u(yy,xx,:,:) = E_temp\u_temp;
			alpha_v(yy,xx,:,:) = E_temp\v_temp;
			iplot = 0; % set iplot to zero after first pass, even if it was 1.
			
%			% weighted least squares... more flexible but 4x as slow.
%			weights = ones(size(Gamma));
%			weights(1:ceil-1) = 1000000;
%			weights = diag(weights); % diagonal matrix of weights
%
%			% estimate modal composition (first 'nmodes' modes only)
%			E_temp = eigenvecs(:,1:nmodes);
%				
%			% u velocity
%			u_temp = squeeze(u(yy,xx,1:ind,:));
%			alpha_u(yy,xx,1:size(E_temp,2),:) = ...
%				inv(E_temp.'*inv(weights)*E_temp)*E_temp.'*inv(weights)*u_temp; 
%				
%			% v velocity
%			v_temp = squeeze(v(yy,xx,1:ind,:));
%			alpha_v(yy,xx,1:size(E_temp,2),:) = ...
%				inv(E_temp.'*inv(weights)*E_temp)*E_temp.'*inv(weights)*v_temp;
%			
%			iplot = 0; % set iplot to zero after first pass, even if it was 1.
		else  % if shallow, set alpha's to nan
			alpha_u(yy,xx,:,:) = nan;
			alpha_v(yy,xx,:,:) = nan;
		end
	end
end

% some statistics for alpha
% root mean (space and time) square for each mode, normalized
alpha_uRMS = zeros(nmodes,1);
alpha_vRMS = zeros(nmodes,1);
for ii = 1:nmodes
	alpha_utemp = squeeze( alpha_u(:,:,ii,:) );
	alpha_uRMS(ii) = sqrt( nanmean( alpha_utemp(:).^2 ) );
	alpha_vtemp = squeeze( alpha_v(:,:,ii,:) );
	alpha_vRMS(ii) = sqrt( nanmean( alpha_vtemp(:).^2 ) );
end
alpha_uRMS = alpha_uRMS./max(alpha_uRMS);
alpha_uRMS
alpha_vRMS = alpha_vRMS./max(alpha_vRMS);
alpha_vRMS

% reconstructed velocities (any time, last spatial point)
tt = 1;
fno = fno+1;
figure(fno);clf;
subplot(121)
plot(squeeze(u(yy,xx,:,tt)),rC)
title('model')
xlabel('u (m/s)')
ylabel('depth (m)')
ax = axis;
hold on
plot([ax(1) ax(2)],[rC_temp(ceil) rC_temp(ceil)],'k:')
hold off
subplot(122)
plot(eigenvecs(:,1:nmodes)*squeeze(alpha_u(yy,xx,:,tt)),rC_temp)
title('reconstructed')
xlabel('u (m/s)')
axis(ax);
hold on
plot([ax(1) ax(2)],[rC_temp(ceil) rC_temp(ceil)],'k:')
hold off

% variance captured by first few modes... at last point in region
a0VAR = var(squeeze(alpha_u(yy,xx,1,:))*eigenvecs(:,1)'); % variance in BT mode, assuming no coupling
a1VAR = var(squeeze(alpha_u(yy,xx,2,:))*eigenvecs(:,2)'); % variance in BC1 mode, assuming no coupling
aVAR = var(squeeze(alpha_u(yy,xx,1:2,:))'*eigenvecs(:,1:2)'); % variance in BT-BC3 modes combined
uVAR = var(squeeze(u(yy,xx,:,:)),0,2); % variance in zonal velocity at this location
fno = fno+1;
figure(fno);clf;
plot(uVAR,rC,'k')
hold on
plot(a0VAR,rC_temp,'b:')
plot(a1VAR,rC_temp,'b--')
plot(aVAR,rC_temp,'b')
hold off
xlabel('variance (m^2/s^2)')
ylabel('depth (m)')
legend('model u','BT','BC1','BT and BC1-3','location','southeast')
title('zonal velocity variance')

% hovmoller diagrams for modes
YY = 1;%floor(length(lat/2));
V = -.1:.01:.1;
fno = fno+1;
figure(fno);clf;
ax(1)=subplot(121);
[c,h] = contourf(lon,time,squeeze(alpha_u(YY,:,1,:))',20);set(h,'linestyle','none');
colormap(jet(20))
colorbar('location','southoutside')
caxis([-.04 .04])
datetick('y',10,'keeplimits','keepticks')
xlabel('longitude')
% ylabel('year')
title('BT mode coefficient  \alpha_{u,0}')
ax(2)=subplot(122);
[c,h] = contourf(lon,time,squeeze(alpha_u(YY,:,2,:))',20);set(h,'linestyle','none');
colormap(jet(20))
colorbar('location','southoutside')
caxis([-.07 .07])
datetick('y',10,'keeplimits','keepticks')
xlabel('longitude')
% ylabel('year')
title('BC mode coefficient  \alpha_{u,1}')
% % make joint colorbar
% h=colorbar('location','southoutside');
% pos=get(ax(2),'position');
% set(h, 'Position', [.13 .11 pos(1)+pos(3)-.13 .0581])
% % reset plot position
% for i=1:2
% pos=get(ax(i), 'Position');
% set(ax(i), 'Position', [pos(1) 2.5*pos(2) pos(3) .8*pos(4)]);
% end

% modal coupling... should make a map of coherence!
YY = 1;
XX = 1;
[r,p] = corrcoef(squeeze(alpha_u(YY,XX,1,:)),squeeze(alpha_u(YY,XX,2,:)))

return

%%

% modal coherence in wavenumber
NW=8;qbias=1;confn=0;qplot=1; 
num = 0;
mean_coherence = zeros(1,floor(size(alpha_u,2)/2));
mean_phase = zeros(1,floor(size(alpha_u,2)/2));
for yy=1%:Nlat 
	yy
	for tt=1:Nreps
		if ~isnan(alpha_u(yy,1,1,tt))
 		[s,coherence,ph,ci,phi] = cmtm(alpha_u(yy,:,1,tt),alpha_u(yy,:,2,tt),dx*distX,NW,qbias,confn,qplot);
 		mean_coherence = mean_coherence + coherence;
		mean_phase = mean_phase + ph;
		num = num+1;
		end
	end
end

% fno = fno+1;
figure(fno);clf;
subplot(2,1,1,'replace')
plot(s,mean_coherence/num,'k')
hold on
plot([min(s) max(s)],[.4 .4],'k--')
hold off
axis tight
ylim([0 1])
xlabel('wavenumber (cpk)')
ylabel('coherence')

subplot(2,1,2,'replace')
plot(s,mean_phase/num,'k')
hold on
plot([min(s) max(s)],[.4 .4],'k--')
hold off
axis tight
ylim([-180 180])
xlabel('wavenumber (cpk)')
ylabel('phase')


%%

% modal coherence in frequency... this is slow.
fno = fno+1;
figure(fno);clf;
NW=8;qbias=1;confn=0;qplot=0; 
num = 0;
dn = 1;
cBIG = nan(length(lat(1:dn:Nlat)),length(lon(1:dn:Nlon)),96);
for yy=1:dn:Nlat 
	yy
	for xx=1:dn:Nlon
		if ~isnan(alpha_u(yy,xx,1,1))
 		[s,cBIG((yy+dn-1)/dn,(xx+dn-1)/dn,:),ph,ci,phi] = cmtm(alpha_u(yy,xx,1,:),alpha_u(yy,xx,2,:),dt,NW);
%		[cBIG((yy+dn-1)/dn,(xx+dn-1)/dn,:),s] = mscohere(squeeze(alpha_u(yy,xx,1,:)),squeeze(alpha_u(yy,xx,2,:)),hamming(round(size(alpha_u,4)/8)),.5,2^nextpow2(size(alpha_u,4)),1/(dt*86400));
 		num = num+1;
% 		drawnow; pause(.1)
		end
	end
end

c = squeeze(sum(sum(cBIG(~isnan(cBIG(:,1,1)),~isnan(cBIG(1,:,1)),:),1),2))/num;
fno = fno+1;
figure(fno);clf;
plot(s,c);
axis([s(1) s(end) 0 1])
hold on
plot(s,ci,'--');
hold off

c_map = mean(cBIG,3);
fno = fno+1;
figure(fno);clf;
hold on
m_proj('Miller Cylindrical','lon',[0 360],'lat',[-66 66]);
%m_proj('Miller Cylindrical','lon',[lon(1)-20 lon(end)+20],'lat',[lat(1) lat(end)]);
[c,h] = m_contourf(lon(1:dn:Nlon),lat(1:dn:Nlat),c_map,20);set(h,'linestyle','none');colorbar
%caxis([])
ca = caxis;
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('averaged coherence map')
% add bathymetry
[cs,h] = m_etopo2('contour',[-6000:2000:0],'edgecolor',[.7 .7 .7]);
%clabel(cs,h,v,'fontsize',6);
caxis(ca);
m_coast('patch','k','edgecolor','k');
m_grid('box','on','tickdir','in');
hold off

%save([basedir 'data/ecco2/cs510_runs/cube84/saved/CoherenceMap_' num2str(latvec(1)) '_' num2str(latvec(2)) '_lon_' num2str(lonvec(1)) '_' num2str(lonvec(2))],'c_map','lat','lon','SaveDate');

%% global coherence map from saved data

% datadir = [basedir 'data/ecco2/cs510_runs/cube84/saved/'];
% 
% % get filenames
% FileNames = dir([datadir 'CoherenceMap_*.mat']);
% % convert to cell and reshape so that each column is one day
% FileNames = {FileNames.name};
% 
% fno = fno+1;
% figure(fno);clf;
% hold on
% ca = [0 1];
% for ii=1:length(FileNames)
%     load([datadir FileNames{ii}]);
%     % plot
%     m_proj('Miller Cylindrical','lon',[0 360],'lat',[-66 66]);
%     [c,h] = m_contourf(lon,lat,c_map,20);set(h,'linestyle','none');colorbar
%     caxis(ca);
% end
% set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
% title('averaged coherence map')
% % add bathymetry
% [cs,h] = m_etopo2('contour',[-6000:2000:0],'edgecolor',[.7 .7 .7]);
% m_coast('patch','k','edgecolor','k');
% m_grid('box','on','tickdir','in');
% hold off
% 
% % should add solid contour at level of no significance

