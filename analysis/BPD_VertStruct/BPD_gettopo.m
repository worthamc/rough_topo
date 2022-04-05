function [h,y] = BPD_gettopo(lat,lon)

% part of the BPD solving set of functions.  Can be called by BPD_map.m or
% other similar function, and returns a particular realization of a
% topography profile with a spectrum similar to the topograhy spectrum near
% (lat,lon) from the ECCO2 model topography ("Dep_ll")

% C. Wortham 21 Aug. 2020

global lon_qd lat_qd topo_x topo_y Dep_ll


 %% depth and y grids

Dep_temp = Dep_ll(abs(lon_qd-lon)<=topo_x , abs(lat_qd-lat)<=topo_y);

dy = 1000*sw_dist([lat,lat],[lon,lon+.25],'km'); % distance (in m) corresponding to 1/4 degree
y = 0:(size(Dep_temp,2)-1);
y = y*dy; % y grid in meters

%% average topography profile

L = size(Dep_temp,2);
NFFT = 2^nextpow2(L);

% nan land... or don't
Dep_temp2 = Dep_temp;
% Dep_temp2(Dep_temp==0) = nan;

% detrend each column
Dep_temp2 = detrend(Dep_temp2')';

% compute average spectrum
Dep_spec = fft(Dep_temp2,NFFT,2);
Dep_spec = 2*abs(Dep_spec(:,1:NFFT/2+1));
Dep_spec = mean(Dep_spec,1);

% wavenumber (cpm)
l = 1/(2*dy)*linspace(0,1,NFFT/2+1);

%     % normalize spectrum to maintain variance
%     A = var(Dep_temp2(:)); % record variance
%     dl = l(2)-l(1);
%     B = sum(Dep_spec*dl); % integral of the normalized spectrum
%     Dep_spec = Dep_spec*A/B; % normalized such that the integral of the spectrum is then record variance.

% kill long wavelength components of topography, satisfying Bobrovich
% and Reznik (1999) eq. 2.7b
Dep_spec( (1./l)>300e3 ) = 0;

% create realization of topography
phi = 2*pi*rand([1,length(Dep_spec)]); % random phase
h = zeros(size(y));
for jj=1:length(Dep_spec)
	h = h + sqrt(Dep_spec(jj)) * cos(2*pi*l(jj)*y + phi(jj));
end

% % maintain std. dev.
% % This maintains the total variance, effectively increasing the amplitude
% % of the small-scale part relative to what was in the original profile.
% h = h * sqrt(var(Dep_temp(:)) / var(h)); 

% modify the amplitude to see the effect (optional)
% h = 2*h;

% DC offset
h = h + nanmean(Dep_temp(:));

%% optional stuff

%     % interpolate to double horizontal grid resolution.
%     rate = 2;
%     y = interp(y,rate);
%     h = interp(h,rate);
% %     Y = linspace(min(y),max(y),rate*length(y)-1);
% %     h = interp1(y,h,Y);
% %     y = Y;

%     % or pick one topography profile
%     Dep_temp2 = Dep_temp;
%     Dep_temp2 = detrend(Dep_temp2')';
%     ind = round(1 + (length(Dep_temp)-1).*rand(1));
%     h = Dep_temp2(:,ind)' + nanmean(Dep_temp(:));



