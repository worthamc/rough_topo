function [eigenvecs_mean,eigenvecs2D,w_r,w_i,eigenvals,indBC] = BPD_sortmodes(S_temp,T_temp,rC_temp,lat,k,y,h,target_vals,num_vals,iplot)
% calls BPD_pde.m to solve the eigenvalue problem for the Rossby wave modes
% over topography and then sorts the resulting modes to find the lowest
% baroclinic and barotropic modes.
% Basically, this just identifies indBC, the index of the lowest baroclinic
% mode.  If BPD_pde didn't find the lowest baroclinic mode, this re-runs
% solving the full eigenvalue problem.

% INPUT:
% S: time mean salinity profile (psu)
% T: time mean in situ demperature profile (deg C)
% depth: of model interfaces (m)
% lat: latitude
% k: zonal wavenumber of large scale wave (cpk)
% y: (regular) horizontal grid (m)
% h: topography profile. same size as y. (m)
% target_vals: frequency of the BTT first baroclinic mode.
% num_vals: number of eigenvectors to solve for.  This just lets me run
%	eigs, which is a faster sparse eigenvalue solver.
% iplot: 1 to draw plots, 0 to skip
%
% OUTPUT:
% eigenvecs_mean: 2D array of large-scale eigenvectors.  This is the mean
%	of eigenvecs2D in across the horizontal dimension.
% eigenvecs2D: 3D array of eigenvectors.  Will have size 
%	[length(depth) length(y) length(w_r)]. The number of eigenvalues is
%	variable because I ignore eigenvectors with negative eigenvalues.
% w_r: real part of eigenvals.  (cycles/day)
% w_i: imaginary part of eigencals.  (cycles/day)
% eigenvals: full eigenvalues. (complex)
% indBC: index of the lowest baroclinic mode
%
% C. Wortham, 21 Aug. 2020
%
% 17 Feb 2021: Solving the eigenvalue problem (BPD_pde.m) with frequency
% omega as the eigenvalue.  Before, had alpha=1/omega as the eigenvalue.
% Now also output the real and imaginary parts of omega, both sorted in
% order of real part.

% set this depending on where script is running
basedir = '/Users/worthamc/Documents/research/';
%basedir = '/data/worthamc/research/';
%addpath(genpath([basedir 'matlab/']),'-end')


% setting up a while loop to allow solving the full eigenvalue problem if
% the sparse solutions don't contain the desired mode
solutions = 'sparse'; % start with sparse eigenvalue solver eigs
while true

if strcmp(solutions,'sparse')
	[eigenvecs_mean,eigenvecs2D,w_r,w_i,eigenvals] = BPD_pde(S_temp,T_temp,rC_temp,lat,k,y,h,num_vals,target_vals);
else
	[eigenvecs_mean,eigenvecs2D,w_r,w_i,eigenvals] = BPD_pde(S_temp,T_temp,rC_temp,lat,k,y,h);
end

% now find the lowest baroclinic mode
% identify the modes with y-dependence (change of sign)
ymode = ~(abs(sum(squeeze(sign(real(eigenvecs2D(1,:,:)))),1))==size(eigenvecs2D(1,:,:),2));
% identify modes with high-mode z-dependence (internal maxima)
zmode = ~(abs(sum(sign(diff(real(eigenvecs_mean),1,1)),1))==(size(eigenvecs_mean,1)-1));
% identify barotropic modes (top amplitude less than twice bottom amplitude)
%     BTmode = (abs(eigenvecs_mean(1,:)./eigenvecs_mean(end,:)) < 2); % weak barotropic mode identification criterion
BTmode = squeeze(any(eigenvecs2D(1,:,:)<eigenvecs2D(end,:,:),2))'; % strong barotropic mode identification criterion
% identify higher baroclinic modes... this is too strict to be useful.
% highmode = squeeze(any(any(diff(eigenvecs2D(1:end-1,:,:),1)>=0,1),2))';
% search modes without y or z dependence
TopAmplitude_search = eigenvecs_mean(1,:);
TopAmplitude_search(ymode | zmode) = nan;
%     BotAmplitude_search = eigenvecs_mean(end,:);
%     BotAmplitude_search(ymode | zmode) = nan;
%     TopVar_search = squeeze(var(eigenvecs2D(1,:,:),0,2));
%     TopVar_search(ymode | zmode | BTmode) = nan;
%     Freq_search = omega_BPD';
%     Freq_search(ymode | zmode | BTmode) = nan;
% find the remaining mode with the largest surface amplitude
[val,indBC] = max(abs(TopAmplitude_search));
%     [val,indBC] = min(abs(BotAmplitude_search));
%     [val,indBC] = min(TopVar_search);
%     [val,indBC] = max(Freq_search);

% if you didn't find a suitable mode, solve the full problem.  But only do
% this once
if strcmp(solutions,'sparse')
	if (isnan(val) || indBC==1)
		solutions = 'full';
		disp('failed to find the lowest barotropic mode... solving full eigenvalue problem')
		continue
	end
end

% if you found the barotropic mode return empty output and the calling
% function will re-run with new topography
if abs(eigenvecs_mean(1,indBC)/eigenvecs_mean(end,indBC)) < 2
	disp('found barotropic mode... repeating')
	eigenvecs_mean = [];
	eigenvecs2D = [];
	w_r = [];
	w_i = [];
	eigenvals = [];
	indBC = [];
	return
elseif w_r(indBC)<target_vals
	disp('found a higher baroclinic mode... repeating')
	eigenvecs_mean = [];
	eigenvecs2D = [];
	w_r = [];
	w_i = [];
	eigenvals = [];
	indBC = [];
	return        
end

break % breaking the while loop to continue to next iteration of for
end

if iplot

% make some basic plots
fno = 80;

% plot mean modes
fno = fno+1;
figure(fno);clf
nn=min(20,length(w_r));
plot(eigenvecs_mean(:,2:nn),rC_temp,'linewidth',1,'color',[.5 .5 .5])
% plot(eigenvecs_mean(:,(~ymode & ~zmode)),rC_temp,'linewidth',1)
hold on
plot(eigenvecs_mean(:,1),rC_temp,'k--','linewidth',2)
plot(eigenvecs_mean(:,indBC),rC_temp,'k','linewidth',2)   
% plot(eigenvecs_mean(:,5),rC_temp,'b','linewidth',2)
% plot(eigenvecs_mean(:,6),rC_temp,'g','linewidth',2)
% plot(eigenvecs_mean(:,7),rC_temp,'m','linewidth',2)
% plot(eigenvecs_mean(:,8),rC_temp,'y','linewidth',2)
plot(zeros(size(rC_temp)),rC_temp,'k:')
hold off
% title(['normalized vertical structures for first ' num2str(nn) ' modes'])
ylabel('depth (m)')
xlabel('P^m_n(z)')
drawnow

% plot bathymetry
fno = fno+1;
figure(fno);clf
% plot(y/1000,Dep_temp)
hold on
plot(y/1000,h,'linewidth',2)
hold off
ylabel('bathymetry (m)')
xlabel('meridional distance (km)')
drawnow

% plot mode vertical structures
fno = fno+1;
figure(fno);clf
mode_no = [1;indBC;5;9]; % choose four modes to plot    
% mode_no = [1;2;3;4]; % choose four modes to plot    
mode_id = {'P^0_0','P^0_1','P^2_0','P^2_1'};
ca = [-2.5 2.5];
load([basedir 'matlab/CW/colormap/anomaly.mat'],'anomaly') % load my anomaly colormap
for jj=1:length(mode_no)
	ax(jj) = subplot(2,2,jj);
   % plot mode
	data = eigenvecs2D(:,:,mode_no(jj)); data(data>max(ca))=max(ca);data(data<min(ca))=min(ca);
	contourf(y/1000,rC_temp/1000,data,20,'linecolor','none');
	set(gcf,'colormap',anomaly)
	caxis(ca)
	title(['mode ' mode_id{jj}])
	h1 = line(y/1000,h/1000,'color','k','linewidth',2);
end
% make joint colorbar
hc=colorbar('location','eastoutside');
pos=get(hc,'position');
set(hc, 'Position', [.9 pos(2) pos(3) .8150])
set(get(hc,'Ylabel'),'String','P^m_n(y,z)')
% reset plot position
pos=get(ax(1), 'Position');
set(ax(1), 'Position', [pos(1) pos(2) .9*pos(3) pos(4)]);
pos=get(ax(3), 'Position');
set(ax(3), 'Position', [pos(1) pos(2) .9*pos(3) pos(4)]);
pos=get(ax(2), 'Position');
set(ax(2), 'Position', [.95*pos(1) pos(2) .9*pos(3) pos(4)]);
pos=get(ax(4), 'Position');
set(ax(4), 'Position', [.95*pos(1) pos(2) .9*pos(3) pos(4)]);
% add axes titles
axes(ax(1))
ylabel('depth (km)')
axes(ax(3))
ylabel('depth (km)')
xlabel('meridional distance (km)')
axes(ax(4))
xlabel('meridional distance (km)')
drawnow

% plot mode vertical structures again... just two now.
fno = fno+1;
figure(fno);clf
mode_no = [5;9]; % choose two modes to plot    
mode_id = {'P^0_0','P^0_1'};
ca = [-2.5 2.5];
load([basedir 'matlab/CW/colormap/anomaly.mat'],'anomaly') % load my anomaly colormap
for jj=1:length(mode_no)
ax(jj) = subplot(1,2,jj);
% plot mode
contourf(y/1000,rC_temp/1000,eigenvecs2D(:,:,mode_no(jj)),20,'linecolor','none');
	set(gcf,'colormap',anomaly)
	caxis(ca)
	% title(['mode ' mode_id{jj}])
	h1 = line(y/1000,h/1000,'color','k','linewidth',2);
end
% make joint colorbar
hc=colorbar('location','eastoutside');
pos=get(hc,'position');
set(hc, 'Position', [.9 pos(2) pos(3) .8150])
set(get(hc,'Ylabel'),'String','P^0_1(y,z)')
% reset plot position
pos=get(ax(1), 'Position');
set(ax(1), 'Position', [pos(1) pos(2) .9*pos(3) pos(4)]);
pos=get(ax(2), 'Position');
set(ax(2), 'Position', [.95*pos(1) pos(2) .9*pos(3) pos(4)]);
% add axes titles
axes(ax(1))
ylabel('depth (km)')
xlabel('meridional distance (km)')
axes(ax(2))
xlabel('meridional distance (km)')
drawnow

end






















