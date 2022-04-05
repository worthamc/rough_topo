% C. Wortham, 15 Sept. 2020

% set this depending on where script is running
basedir = '/Users/worthamc/Documents/research/';
%basedir = '/data/worthamc/research/';
%addpath(genpath([basedir 'matlab/']),'-end')

fno = 20;

load([basedir 'scripts/BPD_VertStruct/BPD_phasespeed_global.mat'])
% load([basedir 'scripts/BPD_VertStruct/BPD_phasespeed_2xTOPO.mat'])
% load('~/Desktop/BPD_phasespeed_global.mat')
load([basedir 'matlab/CW/colormap/standard_cmap.mat'],'standard_cmap')
load([basedir 'matlab/CW/colormap/speed_cmap.mat'],'speed_cmap')
load([basedir 'matlab/CW/colormap/anomaly.mat'],'anomaly')

% make lat and lon vectors for easier handling
lat = lat_map(:,1) -dLatLon/2;
lon = lon_map(1,:)' -dLatLon/2;


% load the SpecChar (AVISO spectral characteristics) stuff
load([basedir 'scripts/spectrum/AVISO/msla_merged/SpecChar_global.mat'])
lat_aviso = [SpecChar_global.latitude];
lat_aviso = reshape(lat_aviso,2,[]);
lat_aviso = mean(lat_aviso,1);
lon_aviso = [SpecChar_global.longitude];
lon_aviso = reshape(lon_aviso,2,[]);
lon_aviso = mean(lon_aviso,1);
% phase speed
phasespeed = [SpecChar_global.phasespeed];
phasespeed = [phasespeed.val];
phasespeed = reshape(phasespeed,3,[]);
phasespeed(phasespeed==inf) = nan;
% gridding the data
lat_g = -60:5:60;
lon_g = 0:5:355;
phasespeed_g = nan(length(lon_g),length(lat_g));
for ii=1:length(lat_aviso)
	yy = find(lat_g==lat_aviso(ii));
	xx = find(lon_g==lon_aviso(ii));
	phasespeed_g(xx,yy) = phasespeed(1,ii);
end


% load phase speeds digitized from Chelton and Schlax (1996) Fig. 5
load([basedir 'scripts/BPD_VertStruct/CheltonSchlax1996_fig5.mat'])


%% phase speed plots

% compare to Fu(2009) Fig. 1
fno = fno+1;
figure(fno);clf
V = [-10 0];
data = BPD_phasespeed;
data(data<min(V)) = min(V);
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
m_contourf([lon' lon'+360],lat,[data data],20,'linecolor','none');
caxis(V)
plot_axes = gca;
cbar_axes = colorbar;
set(get(cbar_axes,'Ylabel'),'String','phase speed (km/day)')
c_map = colormap(jet(20));
% c_map = colormap(standard_cmap);
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('with topography')
cbar_triangles(cbar_axes,plot_axes,c_map,0,1,0,0)

% BTT phasespeed
BTT_phasespeed = BTT_omega./wavenumber;

fno = fno+1;
figure(fno);clf
V = [-10 0];
data = BTT_phasespeed;
data(data<min(V)) = min(V);
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
m_contourf([lon' lon'+360],lat,[data data],20,'linecolor','none');
caxis(V)
plot_axes = gca;
cbar_axes = colorbar;
set(get(cbar_axes,'Ylabel'),'String','phase speed (km/day)')
c_map = colormap(jet(20));
% c_map = colormap(standard_cmap);
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('flat bottom')
cbar_triangles(cbar_axes,plot_axes,c_map,0,1,0,0)

% ratio
ratio = BPD_phasespeed./BTT_phasespeed;

fno = fno+1;
figure(fno);clf
V = [1 2.5];
data = ratio;
data(data>max(V)) = max(V);
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
m_contourf([lon' lon'+360],lat,[data data],20,'linecolor','none');
caxis(V)
plot_axes = gca;
cbar_axes = colorbar;
c_map = colormap(jet(20));
% c_map = colormap(standard_cmap);
set(get(cbar_axes,'Ylabel'),'String','c_{topo}/c_{flat}')
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
% title('phase speed ratio')
cbar_triangles(cbar_axes,plot_axes,c_map,1,0,0,0)

% smooth ratio... may show impact of mid atlantic ridge better.
n = 3; % number of neighboring points to smooth over.
ratio_smooth = convn([ratio ratio],ones(n,n)/(n^2),'same');

fno = fno+1;
figure(fno);clf
V = [1 1.05 1.1 1.2 1.4 1.6 1.8 2 2.5 3.5];
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
[c,h] = m_contour([lon' lon'+360],lat,ratio_smooth,V);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
clabel(c,h,'labelspacing',48);
colormap(jet);
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots

% colormapeditor


%% zonal means
BPD_phasespeed_mean = nanmean(BPD_phasespeed,2);
BTT_phasespeed_mean = nanmean(BTT_phasespeed,2);
% mean_ratio = BPD_phasespeed_mean./BTT_phasespeed_mean;
mean_ratio = nanmean(BPD_phasespeed./BTT_phasespeed,2);
mean_ratio(abs(lat)<8) = nan;

% compare with Chelton et al (2007) Fig. 4
% add my zonal mean phase speeds from SpecChar_global.mat to this figure and the next.
fno = fno+1;
figure(fno);clf
plot(lat,-BPD_phasespeed_mean,'r',lat,-BTT_phasespeed_mean,'b','linewidth',2)
ax=axis; ax(1)=-60; ax(2)=60; ax(4)=30; axis(ax)
xlabel('latitude')
ylabel('westward phase speed (km/day)')
% title('zonal mean phase speed')
% add points for AVISO observations
hold on
% plot(lat_g(lat_g>0),-nanmean(phasespeed_g(:,lat_g>0),1),'k:');
% plot(lat_g(lat_g<0),-nanmean(phasespeed_g(:,lat_g<0),1),'k:');
scatter(lat_aviso((lon_aviso>100)&(lon_aviso<180)),-phasespeed(1,(lon_aviso>100)&(lon_aviso<180)),'k+');
% plot(lat_cs1996,speed_cs1996,'ko');
hold off
% % add some phase speeds from T/P.  See SpecCharPlots.m for some details.
% theta = [10 14 20 24 30 34 40 44 50];
% % 150-190 E (west N Pac)
% gamma1 = [18 12 8 5 4 2.5 1.6 NaN .9];
% % 180-220 E (east N Pac)
% gamma2 = [15 10 7 5 3.5 2.5 1.5 NaN 1.1];
% % 300-340 E (N Atl)
% gamma3 = [16 9 5.5 4 3 2.5 NaN NaN NaN];
% % 155-195 E (west S Pac)
% gamma4 = [20 12 8 5.5 4.5 3 NaN 1.6 NaN];
% % 324-360 E (west S Atl)
% gamma5 = [16 10 8 5 4 3 NaN 1.5 NaN];
% % 240-280 E (east S Pac)
% gamma6 = [14 9 4.5 4 3 2 1.5 1.3 NaN];
% % 55-95 E (Indian)
% gamma7 = [18 12 8.5 6.5 4 3 NaN 1.4 NaN];
% hold on
% plot(theta,gamma1,'r.',...
%     theta,gamma2,'r.',...
%     theta,gamma3,'r.',...
%     -theta,gamma4,'r.',...
%     -theta,gamma5,'r.',...
%     -theta,gamma6,'r.',...
%     -theta,gamma7,'r.')
% hold off
legend('c_{topo}','c_{flat}','AVISO','location','northwest');legend('boxoff')

% compare with Chelton and Schlax (1996) Fig. 5
fno = fno+1;
figure(fno);clf
plot(lat,mean_ratio)
ax=axis; ax(1)=-60; ax(2)=60; ax(3)=0; axis(ax)
xlabel('latitude')
ylabel('c/c_{flat}')
% title({'ratio of observed and theoretical phase speeds', 'to flat bottom theory'})
% title('zonal mean ratio of BPD:BTT phase speed')
% add points for AVISO observations
% BTT_interp = interp1(lat(~isnan(BTT_phasespeed_mean)),BTT_phasespeed_mean(~isnan(BTT_phasespeed_mean)),lat_aviso);
% BTT_interp_aviso = interp2(lon,lat,BTT_phasespeed,lon_g,lat_g')';
% mean_ratio_aviso = phasespeed_g./BTT_interp_aviso;
BTT_interp_aviso = interp2(lon,lat,BTT_phasespeed,lon_aviso,lat_aviso);
mean_ratio_aviso = phasespeed(1,:)./BTT_interp_aviso;
BTT_interp_cs1996 = interp1(lat(~isnan(BTT_phasespeed_mean)),BTT_phasespeed_mean(~isnan(BTT_phasespeed_mean)),lat_cs1996);
hold on
% plot(lat_g(lat_g>0),nanmean(mean_ratio_aviso(:,lat_g>0),1),'k');
% plot(lat_g(lat_g<0),nanmean(mean_ratio_aviso(:,lat_g<0),1),'k');
scatter(lat_aviso((lon_aviso>100)&(lon_aviso<180)),mean_ratio_aviso((lon_aviso>100)&(lon_aviso<180)),'k+');
plot(lat_cs1996,-speed_cs1996./BTT_interp_cs1996','ko')
plot(lat,ones(size(lat)),'k:')
hold off
% % add points for T/P observations
% BTT_interp = interp1(lat(~isnan(BTT_phasespeed_mean)),BTT_phasespeed_mean(~isnan(BTT_phasespeed_mean)),theta);
% hold on
% plot(theta,-gamma1./BTT_interp,'r.',...
%     theta,-gamma2./BTT_interp,'r.',...
%     theta,-gamma3./BTT_interp,'r.',...
%     -theta,-gamma4./BTT_interp,'r.',...
%     -theta,-gamma5./BTT_interp,'r.',...
%     -theta,-gamma6./BTT_interp,'r.',...
%     -theta,-gamma7./BTT_interp,'r.')
% hold off
legend('c_{topo}','AVISO','CS1996','location','southwest');legend('boxoff')

% approximate basin crossing time
% just taking a 100 degree wide basin now.  Rough average Pacific width.
dist = nan(length(lat),1);
for ii=1:length(lat)
	dist(ii) = sw_dist([lat(ii) lat(ii)],[0 100],'km');
end
fno = fno+1;
figure(fno);clf
plot(lat,-dist./BTT_phasespeed_mean/365,'k')
% plot(lat,-dist./BTT_phasespeed_mean/365,'k',lat,-dist./BPD_phasespeed_mean/365,'b')
% ax=axis; ax(1)=-60; ax(2)=60; ax(4)=30; axis(ax)
ylim([0 40])
xlim([0 55])
xlabel('latitude')
ylabel('basin crossing time (years)')
hold on
plot(lat,18*ones(size(lat)),'k:') % compare with 18 year duration of altimetry
hold off

% return

%% vertical structure plots

% flip mean eigenvecs so that they are positive at the top
% (they are undetermined up to sign)
flips = sign(real(BPD_mode(:,:,1)));
flips = repmat(flips,[1,1,size(BPD_mode,3)]);
BPD_mode = BPD_mode.*flips;

% find index of bottom point in vertical structure
% then normalize
% then pick out bottom value
rC = -nc_varget([basedir 'data/OCCA_1x1_v2/annual/DDtheta.0406annclim.nc'],'Depth_c');
rC = double(rC);
dr = -diff(rC);
dr = ([dr(1);dr] + [dr;dr(end)])/2;
BPD_mode_bottom = nan(size(lat_map));
for ii=1:numel(lat_map)
    [xx,yy] = ind2sub(size(lat_map),ii);
    if ~isnan(BPD_mode(xx,yy,1))
        ind = find(~isnan(BPD_mode(xx,yy,:)),1,'last'); % find bottom index
        n = -1/rC(ind) * sum( squeeze(BPD_mode(xx,yy,1:ind).*BPD_mode(xx,yy,1:ind)).*dr(1:ind) ); % compute vertical integral
        BPD_mode(xx,yy,:) = BPD_mode(xx,yy,:)/sqrt(n); % normalize
        BPD_mode_bottom(xx,yy) = BPD_mode(xx,yy,ind); % pick out bottom value
    end
end

fno = fno+1;
figure(fno);clf
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
m_contourf([lon' lon'+360],lat,[BPD_mode(:,:,1) BPD_mode(:,:,1)],20,'linecolor','none');
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('amplitude of mode at surface')
colorbar;
% colormap(jet(20))
colormap(standard_cmap)

fno = fno+1;
figure(fno);clf
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
m_contourf([lon' lon'+360],lat,[BPD_mode_bottom(:,:) BPD_mode_bottom(:,:)],20,'linecolor','none');
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('amplitude of mode at bottom')
colorbar;
% colormap(jet(20))
colormap(standard_cmap)
caxis([-.6 .6])

fno = fno+1;
figure(fno);clf
data = real(BPD_mode(:,:,1))./real(BPD_mode_bottom(:,:));
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
% m_contourf([lon' lon'+360],lat,[data data],20,'linecolor','none');
m_pcolor([lon' lon'+360],lat,[data data]);shading flat
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('ratio of top to bottom amplitude')
colorbar;
% colormap(jet(20))
colormap(standard_cmap)
caxis([-10 10])
% compare this to time-mean of [v(top)/v(bottom)] from ECCO2.  That might
% be a more interesting quantity than just the modal coherence.  Of course,
% the modal coherence is important for showing that the ratio
% [v(top)/v(bottom)] has a chance of being somewhat constant.

% decompose into BTT modes.
alpha0 = nan(size(lat_map));
alpha1 = nan(size(lat_map));
for ii=1:numel(lat_map)
    [xx,yy] = ind2sub(size(lat_map),ii);
    if ~isnan(BPD_mode(xx,yy,1))
        y = squeeze(BPD_mode(xx,yy,~isnan(BPD_mode(xx,yy,:)))); % the BPD mode at this point
        A = squeeze(BTT_mode1(xx,yy,~isnan(BPD_mode(xx,yy,:)))); % BTT BC1 mode at this point
        A = [ones(size(A)) A]; % combine with BTT BT mode
        alpha = A\y; % least squares solution for BTT mode decomposition
        alpha0(xx,yy) = alpha(1);
        alpha1(xx,yy) = alpha(2);        
    end
end

fno = fno+1;
figure(fno);clf
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
% m_contourf([lon' lon'+360],lat,[alpha0 alpha0],20,'linecolor','none');
m_pcolor([lon' lon'+360],lat,[abs(alpha0) abs(alpha0)]);shading flat
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('BTT BT mode coefficient')
colorbar;
% colormap(jet(20))
colormap(standard_cmap)
% caxis([0 100])

fno = fno+1;
figure(fno);clf
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
% m_contourf([lon' lon'+360],lat,[alpha1 alpha1],20,'linecolor','none');
m_pcolor([lon' lon'+360],lat,[abs(alpha1) abs(alpha1)]);shading flat
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('BTT BC1 mode coefficient')
colorbar;
% colormap(jet(20))
colormap(standard_cmap)
% caxis([0 100])

fno = fno+1;
figure(fno);clf
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
% m_contourf([lon' lon'+360],lat,[alpha1 alpha1],20,'linecolor','none');
m_pcolor([lon' lon'+360],lat,[abs(alpha1./alpha0) abs(alpha1./alpha0)]);shading flat
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('ratio of BTT BC1 mode to BTT BT mode coefficient')
colorbar;
% colormap(jet(20))
colormap(standard_cmap)
% caxis([0 100])


%% topography and roughness map
% going for something similar to Gille et al (2000).
% basically RMS of high-pass filtered bathymetry.

% load 1/4 deg lat-lon grid from ECCO2
load([basedir 'scripts/spectrum/ecco2/cs510_runs/cube84/LatLon_qd.mat'])

% load model bathymetry from ECCO2
load([basedir 'data/ecco2/cs510_runs/cube84/grid/ll4/Dep_ll.mat']);
Dep_ll = circshift(Dep_ll,[180*4 0 0]); % shift longitude to run 0--360

fno = fno+1;
figure(fno);clf
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
m_contourf([lon_qd lon_qd+360],lat_qd,[Dep_ll' Dep_ll'],20,'linecolor','none');
% m_pcolor([lon_qd lon_qd+360],lat_qd,Dep_ll');shading flat
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
% title('depth (m)')
cbar_axes = colorbar;
set(get(cbar_axes,'Ylabel'),'String','depth (m)')
colormap(jet(20))
% colormap(standard_cmap)

% make roughness map, generally following Gille et al (2000)
rough = Dep_ll;
% setup filter frequencies
fN = 4/2; % Nyquist frequency... equivalent to 1/2 degree
fcutoff = 1/1.5; % frequency cutoff... equivalent to 1/fcutoff degrees
% high-pass filter
HPfilter = fdesign.highpass('Fst,Fp,Ast,Ap',.7*fcutoff/fN,1.1*fcutoff/fN,20,1);
HPdesign = design(HPfilter);
rough = filtfilt(HPdesign.Numerator,1,rough);
rough = filtfilt(HPdesign.Numerator,1,rough')';
% square
rough = rough.^2;
% low-pass filter
LPfilter = fdesign.lowpass('Fp,Fst,Ap,Ast',1.1*fcutoff/fN,1.3*fcutoff/fN,1,20);
LPdesign = design(LPfilter);
rough = filtfilt(LPdesign.Numerator,1,rough);
rough = filtfilt(LPdesign.Numerator,1,rough')';
% square root
rough = sqrt(abs(rough));

% another way


% plot
fno = fno+1;
figure(fno);clf
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
m_contourf([lon_qd lon_qd+360],lat_qd,[rough' rough'],20,'linecolor','none');
% m_pcolor([lon_qd lon_qd+360],lat_qd,rough');shading flat
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
% title(['bottom roughness at sacles below ' num2str(1/fcutoff) '^\circ (m)'])
cbar_axes = colorbar;
set(get(cbar_axes,'Ylabel'),'String',['bottom roughness at sacles below ' num2str(1/fcutoff) '^\circ (m)'])
colormap(jet(20))
% colormap(standard_cmap)
% caxis([0 500])

fno = fno+1;
figure(fno);clf
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
m_contourf([lon' lon'+360],lat,[h_rms h_rms],20,'linecolor','none');
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
% title('h_{rms} (m)')
colormap(jet(20))
% colormap(standard_cmap)
cbar_axes = colorbar;
set(get(cbar_axes,'Ylabel'),'String','h_{rms} (m)')
% hold on
% m_contourf(lon_map,lat_map,h_rms,[300],'linecolor','k');
% hold off
% caxis([0 500])

fno = fno+1;
figure(fno);clf
plot(lat,nanmean(h_rms,2))
xlabel('latitude')
ylabel('zonal mean h_{rms} (m)')
% axis tight
xlim([-60 60])

% return

%% stratification section

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

S_mean = squeeze(nanmean(S_mean,1));
T_mean = squeeze(nanmean(T_mean,1));

pres = nan(length(lat_hydro),length(rC));
for ii=1:length(lat_hydro)
    pres(ii,:) = sw_pres(-rC,lat_hydro(ii));
end

rho_mean = sw_dens(S_mean,T_mean,pres); % since I'm using potential temp, sw_dens will give potential density

fno = fno+1;
figure(fno);clf
rho_mean_temp = rho_mean-1000;
rho_mean_temp(rho_mean_temp<23) = 23;
rho_mean_temp(rho_mean_temp>30) = 30;
h=subplot(3,1,1);
p = get(h, 'pos');
p(2) = p(2)-.05; p(4) = p(4)+.05;
set(h,'pos',p)
contourf(lat_hydro,rC(rC>=-1100)+5,rho_mean_temp(:,rC>=-1100)',20,'linecolor','none');
caxis([23 30])
set(gca,'xaxislocation','top')
ylabel('depth (m)')
subplot(3,1,[2 3])
contourf(lat_hydro,rC+5,rho_mean_temp',20,'linecolor','none');
caxis([23 30])
plot_axes = gca;
ylabel('depth (m)')
xlabel('latitude')
cbar_axes = colorbar('location','southoutside');
c_map = colormap(summer(20));
set(get(cbar_axes,'xlabel'),'String','\sigma_\theta (kg/m^3)');
drawnow
cbar_triangles(cbar_axes,plot_axes,c_map,0,0,1,1)


%% now plotting stratification section
T_insitu_mean = sw_temp(S_mean,T_mean,pres,0); % converting potential temp to in situ for sw_bfrq
[bfrq,vort,p_ave] = sw_bfrq(S_mean',T_insitu_mean',pres');
p_mean = mean(p_ave,2); % just keeping one column... doesn't change much over latitude

% fno = fno+1;
figure(fno);clf
data = log10(sqrt(bfrq));
% ca = [0 .02];
ca = [-4 -2];
h=subplot(3,1,1);
p = get(h, 'pos');
p(2) = p(2)-.07; p(4) = p(4)+.07;
set(h,'pos',p)
contourf(lat_hydro,-p_mean(p_mean<=1000),data(p_mean<=1000,:),20,'linecolor','none');
ylim([-1000 0])
% ca = caxis;
caxis(ca)
set(gca,'xaxislocation','top')
ylabel('depth (m)')
subplot(3,1,[2 3])
contourf(lat_hydro,-p_mean(p_mean>=1000),data(p_mean>=1000,:),20,'linecolor','none');
caxis(ca)
plot_axes = gca;
ylabel('depth (m)')
xlabel('latitude')
cbar_axes = colorbar('location','southoutside');
c_map = colormap(summer(20));
logcolorbar(cbar_axes)
set(get(cbar_axes,'xlabel'),'String','N (s^{-1})');
drawnow
% cbar_triangles(cbar_axes,plot_axes,c_map,0,0,1,1)

return

%% near-surface stratification map

global lon_hydro lat_hydro hydro_x hydro_y S_mean T_mean rC

% lat-lon grid for map... lower left corners
dLatLon = 2; 
lat_map = -80:dLatLon:78;
lon_map = 0:dLatLon:358;
[lon_map,lat_map] = meshgrid(lon_map,lat_map);

% averaging area for hydrography.  This is the 'radius' of the box in deg.
hydro_x = 2.5;
hydro_y = 2.5;

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
N2_map = nan([size(lat_map),50]); % save near-surface buoyancy frequency

tic
for ii=1:length(lat_map(:)) % using linear indexing
% for ii=find(lon_map>=270.1 & lon_map<=360)'
% for ii=find(lon_map==196)'
% for ii=find(lat_map==26 & lon_map==210) 

% if deeper than 1000 m ( rC(29)=-1007m )... and no land in the area... and
% not at the equator.
if ~isnan( T_mean(lon_hydro==lon_map(ii) , lat_hydro==lat_map(ii) , find(rC<-1000,1,'first') ))...
        && ~any(any( isnan(T_mean( abs(lon_hydro-lon_map(ii))<=hydro_x , abs(lat_hydro-lat_map(ii))<=hydro_y , 1 )) ))...
        && abs(lat_map(ii)) >= 3

	% average hydrography
	[S_temp,T_temp,rC_temp] = BPD_gethydro(lat_map(ii),lon_map(ii));
	
	% plot BTT modes, and use that for wavelength
    path(path,[basedir 'scripts/VertStruct']);
    [LdB,eigenvecsB,rhoB,N2cphB] = BTTmodes(S_temp',T_temp',rC_temp',lat_map(ii),0,5); % second to last argument zero or one for plotting

	% save some other output
	[a,b] = ind2sub(size(lat_map),ii);
	N2_map(a,b,1:length(rC_temp)-1) =  N2cphB;
	
end
end
toc

fno = fno+1;
figure(fno);clf
N2_bottom = nan(size(lat_map));
for ii=1:length(lat_map(:))
	[a,b] = ind2sub(size(lat_map),ii);
	ind_bottom = find(isnan(N2_map(a,b,:)),1,'first');
	ind_bottom = max(1,ind_bottom-1);
	N2_bottom(a,b) = N2_map(a,b,ind_bottom);
end
TMcW = 1+ 2*sqrt(N2_bottom) ./ nanmean(sqrt(N2_map(:,:,10:end)),3);
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
m_contourf([lon' lon'+360],lat,[real(TMcW) real(TMcW)],20,'linecolor','none');
plot_axes = gca;
cbar_axes = colorbar;
c_map = colormap(standard_cmap)
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('Tailleux and McWilliams measure of stratification')
% cbar_triangles(cbar_axes,plot_axes,c_map,0,1,0,0)

fno = fno+1;
figure(fno);clf
V = 0.2:.1:2;
bvfrq = 1./nanmean(sqrt(N2_map(:,:,10:end)),3);
m_proj('miller cylindrical','lon',[0 360+30],'lat',[-80 80]); % or robinson
m_contourf([lon' lon'+360],lat,[bvfrq ],20,'linecolor','none');
plot_axes = gca;
cbar_axes = colorbar;
c_map = colormap(standard_cmap)
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('inverse of mean buoyancy frequency below 100m (cph^{-1})')


plot_axes = gca;
cbar_axes = colorbar;
% c_map = colormap(jet(20));
c_map = colormap(speed_cmap);
m_coast('patch',[.6,.6,.6],'edgecolor',[.6,.6,.6]);
m_grid('box','on','tickdir','in','xtick',[0 90 180 270 360],'ytick',[-60 -30 0 30 60]);
set(findobj('tag','m_grid_color'),'facecolor','none'); % necessary for m_countourf plots
title('phase speed with topography (km/day)')
cbar_triangles(cbar_axes,plot_axes,c_map,0,1,0,0)




