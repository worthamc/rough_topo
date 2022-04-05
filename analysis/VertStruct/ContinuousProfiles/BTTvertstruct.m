% function BTTvertstruct
% Solving Sturm-Liouville boundary value problem
% for the eigenfunctions and eigenvalues of the
% "Basic Textboox Theory" for linear Rossby waves,
% using LiScEig package
%
% C. Wortham

%%%%%%%%%%%%%%%%%%%
% options
%%%%%%%%%%%%%%%%%%%
n = 250;                     % number of vertical levels
errtol = 1.e-9;             % error tolerance, [1.e-9,1.e-3]
theBVP = 'BTT_simpleexp';     % name of function containing the
                            % detailed problem (BCs, stratification)
nplot = 4;                  % number of eigenfunctions to plot

%%%%%%%%%%%%%%%%%%%
% computation
%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Users/stallion/Documents/matlab/LiScEig'))
liscstartup
[Gamma,Phi,Phip,z,wint,wipr,dom]=sp(theBVP,n,errtol);

%%%%%%%%%%%%%%%%%%%
% plotting
%%%%%%%%%%%%%%%%%%%
figure
plot(Phi(:,1:nplot),z);
hold on
plot([0 0],[z(1) z(end)],'k--')
hold off
xlabel('Phi_n')
ylabel('depth (m)')
title(theBVP,'interpreter','none');
M=cell(1,nplot);
for ii=1:nplot
    M{ii} = ['\Gamma_' num2str(ii-1) '=' num2str(Gamma(ii))];
end
legend(M,'location','northeastoutside');







