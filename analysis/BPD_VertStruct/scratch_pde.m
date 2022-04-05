% set this depending on where script is running
basedir = '/Users/worthamc/Documents/research/';
%basedir = '/data/worthamc/research/';
%addpath(genpath([basedir 'matlab/']),'-end')

% variables
% k = .001; % (cpk)
lat = 27;

%depth
% depth = linspace(0,-5000,40); % (m)
% or use ecco2 vertical grid
load([basedir 'data/ecco2/cs510_runs/cube84/grid/cs510/dr.txt'])
depth = [ -cumsum(dr)']+dr'/2; % middle of layer depths

% stratification
S = 35*ones(size(depth));
T = 18*exp(-10*depth/depth(end))-1;
%T = 20*ones(size(depth));
% or use these from samples, averaged from ECCO2 near 210E,25N
% load T_sample
% load S_sample
% T = T_temp';
% S = S_temp';
% depth = depth(~isnan(T));

% % truncate hydrography and depth vector at some mean depth.
% T = T(depth>-3200);
% S = S(depth>-3200);
% depth = depth(depth>-3200);

% topography
y = linspace(0,500,20)*10^3; % (m)
load DH; % saved a cherry-picked topography profile.
dh = 1*DH;
% dh = .005*max(y)/length(y) * ones(size(y)); % (m). a simple large scale slope.
% dh = 100*randn(size(y)); % (m)
% dh = [13.4299
%   -24.1497
%    14.3448
%    32.6047
%     9.7779
%    20.6939
%    14.5377
%    -6.0688
%     5.8774
%   -15.7457]';
h = depth(end) + cumsum(dh); % (m)

% or take from some region of ECCO2
% y = linspace(0,500,20)*10^3;
% load([basedir 'scripts/spectrum/ecco2/cs510_runs/cube84/LatLon_qd.mat'])
% load([basedir 'data/ecco2/cs510_runs/cube84/grid/ll4/Dep_ll.mat']);
% Dep_ll = circshift(Dep_ll,[-180*4 0 0]);
% lonind = find(210 == lon_qd);
% latind = find(25 == lat_qd);
% lonind = find(315 == lon_qd);
% latind = find(25 == lat_qd);
% h = Dep_ll(lonind,latind:latind+length(y)-1);
% dh = 1*[0 diff(h)];
% h = depth(end) + cumsum(dh);

% plot BPD solutions, and use deformation radius for wavenumber
path(path,[basedir 'scripts/VertStruct']);
[LdB,eigenvecsB,rhoB,N2cphB] = BTTmodes(S,T,depth,lat,1,5);
k = 1/LdB(2); % (cpk)
k = -k/2

tic
[eigenvecs_mean,eigenvecs2D,alpha,eigenvals] = BPD_pde(S,T,depth,lat,k,y,h);
toc


%% sample diagnostic plots


fno = 3; % initial figure number

fno = fno+1;
figure(fno);clf
plot(y,h)
% hold on
% plot(y,dh,'g')
% hold off
title('topography profile (m)')
xlabel('position (m)')

fno = fno+1;
figure(fno);clf
plot(1./alpha)
title('frequency eigenvalues')

fno = fno+1;
figure(fno);clf
nn=12;
pcolor(y,depth,eigenvecs2D(:,:,nn));shading flat;colorbar
title(['\psi   for mode ' num2str(nn)])

fno = fno+1;
figure(fno);clf;
ColorOrder = get(gcf,'DefaultAxesColorOrder');
ColorOrder = repmat(ColorOrder,[10 1]);
hold on
for nn=1:5
    plot(squeeze(eigenvecs2D(:,:,nn)),depth,':','color',ColorOrder(nn,:))
%     plot(sum(eigenvecs2D(:,1,5*(nn-1)+1:5*(nn-1)+5),3),depth)
end
plot(squeeze(mean(eigenvecs2D(:,:,1:nn),2)),depth,'linewidth',2)
hold off
title(['vertical structure for first ' num2str(nn) ' modes'])

fno = fno+1;
figure(fno);clf
nn=5;
plot(eigenvecs_mean(:,1:nn),depth,'linewidth',2)
hold on
plot(zeros(size(depth)),depth,'k:')
hold off
title(['normalized vertical structures for first ' num2str(nn) ' modes'])
ylabel('depth (m)')
xlabel('\Phi_n')

return

%% analysis picture
% looking at relationship between vertical/horizontal mode number and
% frequency

num=length(alpha);

% count number of sign transitions (zero crossings) in vertical for first
% num modes
z_zeros = nan(1,num);
for nn=1:num
    temp = sign(eigenvecs2D(:,1,nn));
    z_zeros(nn) = sum(abs(diff(temp)))/2;
end

% count number of sign transitions (zero crossings) in horizontal for first
% num modes
y_zeros = nan(1,num);
for nn=1:num
    temp = sign(eigenvecs2D(1,:,nn));
    y_zeros(nn) = sum(abs(diff(temp)))/2;
end

fno = fno+1;
figure(fno);clf
scatter(y_zeros,z_zeros,100,log10(1./alpha(1:num)),'filled');
cbar=colorbar;
logcolorbar(cbar)
xlabel('horizontal mode number (m)')
ylabel('vertical mode number (n)')
title('frequency as function of mode number')
% don't want to contour plot here because the modes are intrinsically
% discrete.  Scatter plot emphasizes this, while contour would obscure the
% fact.

omega1 = 1/alpha(1)


% return

%% dispersion relation


% 1/LdB(2) is the BC1 deformation wavenumber, in radians/km
% note that I'm looking at negative k's.
kmax = -.0065;
Nk = 50;
K = linspace(kmax,kmax/100,Nk); % cpk

% solve BPD_pde at each wavenumber
eigenvecs_mean_disp = nan([size(eigenvecs_mean),Nk]);
eigenvecs2D_disp = nan([size(eigenvecs2D),Nk]);
alpha_disp = nan(length(alpha),Nk);
for kk=1:Nk
    [eigenvecs_mean,eigenvecs2D,alpha] = BPD_pde(S,T,depth,lat,K(kk),y,h);
    eigenvecs_mean_disp(:,:,kk) = eigenvecs_mean;
    eigenvecs2D_disp(:,:,:,kk) = eigenvecs2D;
    alpha_disp(:,kk) = alpha;
end

fno = fno+1;
figure(fno);%clf
hold on
plot(K,2*pi./alpha_disp(1,:))
xlabel('wavenumber (cpk)')
ylabel('frequency (?)')
title('lowest vertical mode dispersion relation')

% add BTT dispersion curves
Omega = 7.292e-5; % rad/s
EarthRadius = 6370e3; % m
beta = @(lat) (2*Omega/EarthRadius)*cosd(lat);
omega_BTT_0 = @(k,lat) -beta(lat) ./ (k.*2*pi/1000);
omega_BTT_1 = @(k,lat) -beta(lat)*(k*2*pi/1000) ./ ( (k*2*pi/1000).^2 + (1./LdB(2)/1000)^2 );

plot(K,omega_BTT_0(K,lat),'b--')
plot(K,omega_BTT_1(K,lat),'b--')
hold off

% frequency output 1/alpha is in cycles/second.  to convert to cycles/day,
% multiply by 86400 seconds/day
