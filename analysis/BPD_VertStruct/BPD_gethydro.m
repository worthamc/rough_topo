function [S_temp,T_temp,rC_temp] = BPD_gethydro(lat,lon)

% part of the BPD solving set of functions.  Can be called by BPD_map.m or
% other similar function, and returns S and T profiles averaged over some
% region near (lat,lon).  Hydrography is based on the OCCA atlas.
%
% assumes input (global) temperature is potential temperature, and 
% outputs in situ temperature.

% C. Wortham 21 Aug. 2020

global lon_hydro lat_hydro hydro_x hydro_y S_mean T_mean rC



%% average hydrography
% pick out all the hydrography in my averaging area and average that above
% the median depth

% pick out the profiles
T_temp = T_mean(abs(lon_hydro-lon)<=hydro_x , abs(lat_hydro-lat)<=hydro_y , :);
S_temp = S_mean(abs(lon_hydro-lon)<=hydro_x , abs(lat_hydro-lat)<=hydro_y , :);

% choose maximum depth
rC_temp = rC(squeeze(any(any(~isnan(T_temp),1),2)));
% and truncate the hydro_temp
T_temp = T_temp(:,:,1:length(rC_temp));
S_temp = S_temp(:,:,1:length(rC_temp));

% find any places with incomplete profiles and nan them... incomplete
% profiles screw up the averaging, creating kinks in the profile
bad_ind = isnan(T_temp(:,:,end));
bad_ind = squeeze(bad_ind);
for jj = 1:numel(bad_ind)
	[a,b] = ind2sub(size(bad_ind),jj);
	if bad_ind(a,b)
	T_temp(a,b,:) = nan; % nan the partial profiles
	S_temp(a,b,:) = nan; % nan the partial profiles
	end
end

% in ECCO2, some salinity profiles have anomalies at the bottom.
% Suspect this is a masking problem.  If the bottom point deviates from
% that above, nan.
% But this isn't quite aggressive enough for ECCO2 hydrography.
% Sometimes the errors persist higher up in the water column.
for jj = 1:numel(S_temp(:,:,1))
	[a,b] = ind2sub(size(S_temp(:,:,1)),jj);
	if abs( S_temp(a,b,end)-S_temp(a,b,end-1) )>0.1
	S_temp(a,b,end) = S_temp(a,b,end-1);
	end
end

% average profiles
T_temp = squeeze(nanmean(nanmean(T_temp,1),2));
S_temp = squeeze(nanmean(nanmean(S_temp,1),2));

% convert to in situ temp
pres = sw_pres(-rC_temp,lat);
T_temp = sw_temp(S_temp,T_temp,pres,0); clear pres;

% % having problems with ecco2 hydrography.  so I'll just replace the last
% % point with a linear extrapolation based on the points above.
% S_temp = interp1(rC_temp(1:end-1),S_temp(1:end-1),rC_temp,'linear','extrap');
% T_temp = interp1(rC_temp(1:end-1),T_temp(1:end-1),rC_temp,'linear','extrap');

%% optional stuff

%     % interpolate to double vertical grid resolution.
%     rate = 2;
%     RC_TEMP = interp(rC_temp,rate);
%     RC_TEMP = RC_TEMP(RC_TEMP>-4500);
%     RC_TEMP = [RC_TEMP(1:end-1)' , RC_TEMP(end):-50:rC_temp(end)]';
%     S_temp = interp1(rC_temp,S_temp,RC_TEMP,'linear','extrap');
%     T_temp = interp1(rC_temp,T_temp,RC_TEMP,'linear','extrap');
%     rC_temp = RC_TEMP;

%     % this was just for testing
%     Dep_temp = min(Dep_temp(:))*ones(size(Dep_temp));
%     S_temp = mean(S_temp)*ones(size(S_temp));
%     T_temp = mean(T_temp)*ones(size(T_temp));