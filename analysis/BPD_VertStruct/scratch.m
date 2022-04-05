% set this depending on where script is running
basedir = '/Users/worthamc/Documents/research/';
%basedir = '/data/worthamc/research/';
%addpath(genpath([basedir 'matlab/']),'-end')

fno = 0; % initial figure number

% variables
k = .001; % (cpk)
l = 0; % (cpk)
lat = 40;

% stratification
depth = 0:-50:-4000; % (m)
S = 35*ones(size(depth));
T = 20*exp(-10*depth/depth(end))-1;
%T = 20*ones(size(depth));

% bottom topography spectrum
lambda_m = linspace(10,100); % (km)
% n=1:100;n=fliplr(n);
% lambda_m = 100./n; % (km)
lambda_m = [-fliplr(lambda_m) lambda_m]; % double, because that's the way I did the theory
lm = 1./lambda_m; % (cpk)

% lm = n.*(1/100); % (cpk)
% lm = [-fliplr(lm) lm]; % double, because that's the way I did the theory
% lambda_m = 1./lm; % (km)

hm = 15*log(abs(lambda_m)/max(lambda_m)*15); % (m)
%hm = 5*log(abs(lambda_m)/max(lambda_m)*1000); % (m)
%hm = 30*ones(size(lm));
% introduce some random phase
hm = hm.*exp(1i*2*pi*rand(size(hm)));

% call LargeScalePsi_z
[psi_bar,Gamma,D,psi_tilde,N2cph] = LargeScalePsi_z(S,T,depth,lat,k,l,lm,hm);

%% Dispersion relation

global G F RHO0 BETA

% choose range of wavenumbers
path(path,[basedir 'scripts/VertStruct']);
[LdB,eigenvecsB,rhoB,N2cphB] = BTTmodes(S,T,depth,lat,1,5);
k_d = 1/LdB(2)/2/pi; % (cpk)

% 1/LdB(2) is the BC1 deformation wavenumber, in radians/km
kmax = 4*k_d;
nk = 50;
nl = 50;
K = linspace(-kmax,kmax,nk); % cpk
%L = linspace(-kmax,kmax,nl); % cpk
L=0; nl=1;
[K,L] = meshgrid(K,L);

vecs_bar = nan(nl,nk,length(depth),length(Gamma));
vecs_tilde = nan(nl,nk,length(depth),length(lm)/2);
gamma = nan(nl,nk,length(Gamma));
for kk=1:nk
    for ll = 1:nl
        [psi_bar,Gamma,D,psi_tilde] = LargeScalePsi_z(S,T,depth,lat,K(kk),L(ll),lm,hm);
        vecs_bar(ll,kk,:,:) = psi_bar;
        vecs_tilde(ll,kk,:,:) = psi_tilde;
        gamma(ll,kk,:) = Gamma;
    end
end

nmode = 2;
Kradm = 2*pi.*K./1000; % (rad/m)
omega = -BETA.*Kradm(1,:) ./ (Kradm(1,:).^2 + gamma(1,:,nmode));

%% figures

fno = 0;

fno = fno+1;
figure(fno)
plot(lambda_m(lambda_m>=0),abs(hm(lambda_m>=0)))
title('topography spectrum')
xlabel('wavelength (km)')
ylabel('amplitude (m)')
% describe as "wavelength spectrum in amplitude"

fno = fno+1;
figure(fno)
plot(lm(lm>=0),abs(hm(lm>=0)))
title('topography spectrum')
xlabel('wavenumber (cpk)')
ylabel('amplitude (m)')
% describe as "spectrum in amplitude"

% reconstruct topography
y = linspace(0,1/min(abs(lm)),2*length(lm));
h = zeros(size(y));
for ii=1:length(lm)
h = h + hm(ii)*exp(1i*2*pi*y*lm(ii));
end
fno = fno+1;
figure(fno)
plot(y,real(h))
title('sample topography')
xlabel('height (m)')
ylabel('distance(km)')

fno = fno+1;
figure(fno)
hold on
plot(K(1,:),omega)
yl = ylim;
ylim([0 yl(2)])
title('dispersion relation')
xlabel('k (cpk)')
ylabel('\omega (cpd)')
hold off

lplot = 1;
[val,kplot] = min(abs(K + k_d/2));

fno = fno+1;
figure(fno)
plot(squeeze(vecs_bar(lplot,kplot,:,1:5)),depth)
title(['large scale vertical structure \Phi at k = ' num2str(K(kplot),2) ' cpk'])
xlabel('amplitude')
ylabel('depth (m)')

fno = fno+1;
figure(fno)
plot(squeeze(vecs_tilde(lplot,kplot,:,:)),depth)
title(['small scale vertical structure \psi at k = ' num2str(K(kplot),2) ' cpk'])
xlabel('amplitude')
ylabel('depth')


