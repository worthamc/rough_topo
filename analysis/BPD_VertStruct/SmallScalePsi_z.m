function [psi_tilde] = SmallScalePsi_z(S,T,k,lm,depth,lat)

% Solving for vertical structure of small horizontal scale streamfunction
% in the "bottom pressure decoupling" scenario, with a large scale wave
% over small scale topography.  Follows Bobrovich and Reznik (1999), though
% modified to apply to an arbitrary stratification.
%
% Equation is in form
% Dzz psi_tilde - lm^2 psi_tilde = 0

% INPUT:
% S: time mean salinity profile (psu)
% T: time mean in situ demperature profile (deg C)
% k: zonal wavenumber of large scale wave (rad/m)
% lm: meridional wavenumber of small scale topography component (rad/m)
% depth: of model interfaces (m)
% lat: latitude
%
% OUTPUT:
% psi_tilde: vertical structure of small scale wave (m^2/s)
%
% C. Wortham, 7 Oct. 2019

% global constants
global G F RHO0

% calculate layer thicknesses
dr = -diff(depth);
dr = ([dr(1) dr] + [dr dr(end)])/2;

% calculate pressure from depth
pres = sw_pres(-depth,lat);

% make sure it's just a vector
dr = squeeze(dr);
N = length(dr);
if (length(S)~=N)
    error('S and dr must be the same length');
end

% make sure that S, T and dr are row vectors
if (size(dr,2)==1)
    dr = dr';
end
if (size(S,2)==1)
    S = S';
end
if (size(T,2)==1)
    T = T';
end

% Dzz term

% create some reuseable matrices
n_minus1 = diag(ones(1,N-1),-1);
n = diag(ones(1,N));
n_plus1 = diag(ones(1,N-1),1);

% reference pressures
pr_left = ( pres([2 2:N]) + pres([1 1:N-1]) )/2;
pr_right = ( pres([2:N N]) + pres([1:N-1 N-1]) )/2;

% denominator, e.g. rho(n) - rho(n-1), for differential operator 
delta_Rho_left = sw_pden(S([2 2:N]),T([2 2:N]),pres([2 2:N]),pr_left) - sw_pden(S([1 1:N-1]),T([1 1:N-1]),pres([1 1:N-1]),pr_left);
delta_Rho_right = sw_pden(S([2:N N]),T([2:N N]),pres([2:N N]),pr_right) - sw_pden(S([1:N-1 N-1]),T([1:N-1 N-1]),pres([1:N-1 N-1]),pr_right);
% deal with negative gradients
if any(delta_Rho_left<=0)
   if delta_Rho_left(end)<=0
       delta_Rho_left(end) = 1e-5; % make bottom positive
   end
   ind = find(delta_Rho_left<=0);
   for ii=1:length(ind)
       delta_Rho_left(ind(ii)) = delta_Rho_left(ind(ii)+1); % work your way up, replacing with next lower value
   end
end
if any(delta_Rho_right<=0)
   if delta_Rho_right(end)<=0
       delta_Rho_right(end) = 1e-5; % make bottom positive
   end
   ind = find(delta_Rho_right<=0);
   for ii=1:length(ind)
       delta_Rho_right(ind(ii)) = delta_Rho_right(ind(ii)+1); % work your way up, replacing with next lower value
   end
end
% make matrices
delta_Rho_left_matrix = repmat(delta_Rho_left',[1 N]);
delta_Rho_right_matrix = repmat(delta_Rho_right',[1 N]);

dr_matrix = repmat(dr',[1 N]);

Dzz = dr_matrix.^-1 .* ( (n_minus1 - n)./delta_Rho_left_matrix - (n - n_plus1)./delta_Rho_right_matrix);
% top boundary condition for Dzz term
Dzz([1 end],:) = zeros(2,N);
Dzz(1,1:2) = 2*[-1 1] / (dr(1) * delta_Rho_right(1));
Dzz(N,[N-1 N]) = 2*[1 -1] / (dr(N) * delta_Rho_left(N));

% lm^2 term
lm2 = lm^2 * eye(N);

% LHS
A = (F^2*RHO0/G) * Dzz - lm2;

% RHS
d = zeros(N,1);
% bottom boundary condition through d term
d(end) = F^3*RHO0 / (k*G*delta_Rho_left(N));

% solve
psi_tilde = A \ d;
