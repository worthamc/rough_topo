function [Ld,eigenvecs,rho,N2cph] = BTTmodes(S,T,depth,lat,iplot,nplot)

% function [Ld,eigenvecs,rho,N2cph] = BTTmodes(S,T,depth,lat,iplot,nplot)
%
% Compute eigenvectors and eigenvalues for the
% Sturm-Liouville problem
% d/dz(f^2/N^2 d/dz(phi)) + Gamma phi = 0
% with BCs
% d/dz(phi) = 0 at z=0,-D
% phi(z) is the vertical structure corresponding
% to "basic textbook theory" Rossby waves 
% linearized about a resting basic state.
%
% Takes S and T profiles as input, computes the
% density profile and buoyancy frequency using 
% a neutral density calculation and then solves 
% discretized eigenvalue problem.
%
% Follows method laid out in Smith (2007, Appendix B)
%
% INPUT:
% S: time mean salinity profile (psu)
% T: time mean in situ demperature profile (deg C)
% depth: of model interfaces (m)
% lat: latitude
% iplot: 1 to plot some output, 0 to skip
% nplot: number of modes to plot, <= length(depth)
%
% OUTPUT:
% Ld: deformation radius (km)
% eigenvecs: associated amplitude profiles
% rho: potential density profile (kg/m^3)
% N2cph: buoyancy frequency squared (cph^2)
% 
% C. Wortham, 24 Nov. 2018

% % sample input
% lat = 25;
% depth = 0:-50:-4000;
% S = 35*ones(size(depth));
% T = 6*(atan(-4*pi*depth/depth(end)+1.5)+3);
% % T = 1 + 30*exp(depth/600);
% % S = linspace(37,35,length(depth));
% % T = linspace(30,2,length(depth));
% iplot = 1;
% nplot = 5;

% hopefully the only NaNs are at the end
depth(isnan(T)) = [];
S(isnan(T)) = [];
T(isnan(T)) = [];

% constants
g = sw_g(lat,0);
f = sw_f(lat);
rho0 = 1000;

% calculate distance between midpoints
dr = -diff(depth);
dr = ([dr(1) dr] + [dr dr(end)])/2;

% calculate pressure from depth
pres = sw_pres(-depth,lat);

% calculate vertial structure function
Lz = (f^2 * rho0 / g) * BTTmatrix(S,T,pres,dr);

% solve problem
[eigenvecs,eigenvals] = eig(Lz);

% transform and sort eigenvalues... note that eigenvectors are not in their original order!!!
Gamma = -diag(eigenvals);
[Gamma, ix] = sort(Gamma);
Ld = 1./real(sqrt(Gamma))/1000; % km
eigenvecs = eigenvecs(:,ix);

% flip eigenvecs so that they are positive at the top
% (they are undetermined up to sign)
flips = sign(eigenvecs(1,:));
flips = repmat(flips,[length(eigenvecs),1]);
eigenvecs = eigenvecs.*flips;

% normalize
% dr_mat = repmat(dr',[1,size(eigenvecs,1)]);
% n = -1/depth(end)*sum(eigenvecs.*eigenvecs.*dr_mat,1);
n = 1/(depth(end))* trapz(depth,eigenvecs.^2,1);
n_mat = repmat(n,[size(eigenvecs,1),1]);
eigenvecs = eigenvecs./sqrt(n_mat);

% compute potential density profile
rho = sw_pden(S,T,pres,0);

% compute N2
N = length(dr);
pr_center = ( pres(2:N) + pres(1:N-1) )/2;
dz = depth(2:N) - depth(1:N-1);
delta_Rho = sw_pden(S(2:N),T(2:N),pres(2:N),pr_center) - sw_pden(S(1:N-1),T(1:N-1),pres(1:N-1),pr_center);
N2 = -(g/rho0) * delta_Rho./dz;
N2cph = (3600/2/pi)^2 * N2; % N2 in cph^2

% plot
if iplot == 1
    figure(99);clf;
    
    subplot(1,4,1:2)
%     plot(eigenvecs(:,1:nplot),depth,'k','linewidth',2);
	set(gca,'LineStyleOrder',{'-','--','-.',':',':x',':o',':^',':+'})
	hold all
	for ii=1:nplot
		plot(eigenvecs(:,ii),depth/1000,'k','linewidth',1.5,'markersize',5);
	end
    plot([0 0],[depth(1) depth(end)]/1000,'k-')
    hold off
%     xlabel('F_n(z)')
    ylabel('depth (km)')
    M=cell(1,nplot);
    for ii=1:nplot
        M{ii} = sprintf('L_d = %0.3g km',Ld(ii));
    end
    legend(M,'location','southwest');
    legend('boxoff')
    axis tight
	yl = ylim;
	box on
%     ylim auto
    
    subplot(143)
    plot(rho-1000,depth/1000,'k')
    xlabel('\sigma_{\theta}  (kg m^{-3})')
    set(gca,'yticklabel',[])
%     a143 = axis;
	ylim(yl)
    
    subplot(144)
    plot(sqrt(N2cph),-sw_dpth(pr_center,lat)/1000,'k')
    xlabel('N  (cph)')
    set(gca,'yticklabel',[])
%     a144 = axis; a144(3) = a143(3);
%     axis(a144)
	ylim(yl)
    
    drawnow
end



