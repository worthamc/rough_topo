function [psi_bar,Gamma,D,psi_tilde,N2] = LargeScalePsi_z(S,T,depth,lat,k,l,lm,hm)

% Solving for vertical structure of large horizontal scale streamfunction
% in the "bottom pressure decoupling" scenario, with a large scale wave
% over small scale topography.  Follows Bobrovich and Reznik (1999), though
% modified to apply to an arbitrary stratification. 
%
% Solves eigenvalue problem for the vertical modes psi_bar and eigenvalues
% Gamma. It turns out to be a Quadratic Eigenvalue Problem (involves square
% of the eigenvalue) of the form
% ( lambda^2 C + lambda B + A ).x = 0
% For more information, see "Templates for the Solution of Algebraic
% Eigenvalue Problems" (Z. Bai et. al, 2000) section 9.2. or "The Quadratic
% Eigenvalue Problem" (F. Tisseur and K. Meerbergen, SIAM, 2001).
%
% Note that the input arguments here are in cpk and cpd, since that's what
% I use in the spectra.  But everything is converted to mks units for this
% calculation and for passing to Dcalc.m and SmallScalePsi_z.m

% INPUT:
% S: time mean salinity profile (psu)
% T: time mean in situ demperature profile (deg C)
% depth: of model interfaces (m)
% lat: latitude
% k: zonal wavenumber of large scale wave (cpk)
% l: meridional wavenumber of large scale wave (cpk)
% lm: meridional wavenumber of small scale topography component (cpk)
% hm: amplitude of each topography wavenumber (m)
%
% OUTPUT:
% psi_bar: vertical structure of large scale waves (m^2/s)
% Gamma: eigenvalues... related to the deformation radii
%
% C. Wortham, 7 Oct. 2019

% change to mks units up front
k = 2*pi*k/1000; % (rad/m)
l = 2*pi*l/1000; % (rad/m)
lm = 2*pi*lm/1000; % (rad/m)
K = sqrt( k^2 + l^2 ); % (rad/m)

% global constants (mks units)
global G F RHO0 BETA
G = sw_g(lat,0);
F = sw_f(lat);
RHO0 = 1000;
BETA = 2*7.292e-5*cosd(lat)/6370e3;

% calculate D
[D,psi_m_tilde] = Dcalc(S,T,depth,lat,k,lm,hm);

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

% matrix A
Dzz = dr_matrix.^-1 .* ( (n_minus1 - n)./delta_Rho_left_matrix - (n - n_plus1)./delta_Rho_right_matrix);
% boundary conditions
Dzz([1 end],:) = zeros(2,N);
Dzz(1,1:2) = 2*[-1 1] / (dr(1) * delta_Rho_right(1));
Dzz(N,[N-1 N]) = 2*[1 -( 1+k*G^2*D*delta_Rho_left(N)^2*K^4 / ( dr(N)*BETA^2*F^3*RHO0^2 ) ) ] / (dr(N) * delta_Rho_left(N));
A = (F^2*RHO0/G) * Dzz;

% matrix B
B = eye(N);
% bottom BC
B(N,N) = 1 + ( 4*k*G*D*K^2*delta_Rho_left(N) ) / ( dr(N)^2*BETA^2*F*RHO0 );

% matrix C
C = zeros(N);
% bottom BC
C(N,N) = ( 2*k*G*D*delta_Rho_left(N) ) / ( dr(N)^2*BETA^2*F*RHO0 );

% if D=0 then C=0, and it's a regular eigenvalue problem.  Otherwise, use
% polyeig.
if D == 0
    [eigenvecs,eigenvals] = eig(A,-B);
else
    [eigenvecs,eigenvals,s] = polyeig(A,B,C);
end

% sort output by eigenvalue
if ~(size(eigenvals,2)==1) % eig returns a matrix of eigenvalues, while polyeig returns a vector.  Fix this.
    eigenvals = diag(eigenvals);
end
eigenvals = eigenvals(eigenvals>=0); % neglect negative eigenvalues... do I also neglect the corresponding eigenvectors?
[Gamma,ind] = sort(eigenvals);
eigenvecs = eigenvecs(:,ind);

% flip eigenvecs so that they are positive at the bottom
% (they are undetermined up to sign)
flips = sign(real(eigenvecs(end,:)));
flips = repmat(flips,[size(eigenvecs,1),1]);eigenvecs = eigenvecs.*flips;
 
% normalize eigenvectors so they integrate to 1 over depth
dr_mat = repmat(dr',[1,size(eigenvecs,2)]);
n = -1/depth(end)*sum(eigenvecs.*conj(eigenvecs).*dr_mat,1);
n_mat = repmat(n,[size(eigenvecs,1),1]);
eigenvecs = eigenvecs./sqrt(n_mat);
% alternative: make maximum to one
% for ii = 1:size(eigenvecs,2)
%     eigenvecs(:,ii) = eigenvecs(:,ii)./max(eigenvecs(:,ii));
% end

psi_bar = eigenvecs;

% calculate buouyancy frequency (needed for psi_tilde)
pr_center = ( pres(2:N) + pres(1:N-1) )/2;
dr = depth(2:N) - depth(1:N-1);
delta_Rho = sw_pden(S(2:N),T(2:N),pres(2:N),pr_center) - sw_pden(S(1:N-1),T(1:N-1),pres(1:N-1),pr_center);
N2 = -(G/RHO0) * delta_Rho./dr; % (rad/s)
N2cph = (3600/2/pi)^2 * N2; % N2 in cph^2

% now go on and calculate the full psi_tilde
% currently just using Gamma(1)... but probably want that associated with
% the lowest mode?  or something.  Didn't realize that these would have
% modes, but it comes in through the frequency.
ymax = 2*pi/min(abs(lm)); % maximum wavelength (m)
y = linspace(0,ymax);
psi_tilde = zeros(length(depth),length(y));
nn = 1; % choose vertical mode
for ii=1:length(lm) % sum over each topography component...
    psi_tilde = psi_tilde + lm(ii)*hm(ii)*psi_m_tilde(:,ii)*exp(1i*lm(ii)*y);
end
psi_tilde = real(1i*k^2*N2(end)*psi_bar(end,nn)./(( -BETA*k/(K^2+Gamma(nn)) )*F^2)*psi_tilde);


