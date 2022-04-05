function [D,psi_m_tilde] = Dcalc(S,T,depth,lat,k,lm,hm)

% Calulates "D", the coupling coefficient between the large scale and small
% scale waves.  D appears in the bottom boundary condition for the large
% scale wave. Solving for vertical structure of large horizontal scale
% streamfunction in the "bottom pressure decoupling" scenario, with a large
% scale wave over small scale topography.  Follows Bobrovich and Reznik
% (1999), though modified to apply to an arbitrary stratification.

% INPUT:
% S: time mean salinity profile (psu)
% T: time mean in situ demperature profile (deg C)
% depth: of model interfaces (m)
% lat: latitude
% k: zonal wavenumber of large scale wave (rad/m)
% omega: frequency of the large scale wave (rad/s)
% lm: meridional wavenumber of small scale topography component (rad/m)
% hm: amplitude of each topography wavenumber (m)
%
% OUTPUT:
% D (m^2/s)
%
% C. Wortham, 7 Oct. 2019

% global constants
global G F RHO0

% calculate psi_tilde(-H) for each wavenumber
psi_tilde_bot = nan(size(lm));
psi_m_tilde = nan(length(depth),length(lm));
for ii=1:length(psi_tilde_bot)
    [psi_tilde] = SmallScalePsi_z(S,T,k,lm(ii),depth,lat);
    psi_tilde_bot(ii) = psi_tilde(end);
    psi_m_tilde(:,ii) = psi_tilde;
end

% calculate D
D = sum( lm.^2 .* abs(hm).^2 .* abs(psi_tilde_bot) );

