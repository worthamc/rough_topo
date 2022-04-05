function [eigenvecs_mean,eigenvecs2D,w_r,w_i,eigenvals] = BPD_pde(S,T,depth,lat,k,y,h,varargin)

% Solve the 2-d (y-z) pde for the vertical structure of Rossby waves over rough
% topography.  See notes on topography and Bobrovich and Reznik (1999) for
% details.  The problem is turns into an eigenvalue pde problem for the
% vertical structure and frequency.  Discretizing and transforming from
% array to linear indices allows me to write the equation as a matrix
% eigenvalue problem of the form
% omega*(A+B+C)*psi + D*psi = 0
% where A represents the meridional derivatives, B represents the vertical
% derivatives, C comes from zonal derivative, and D includes the
% beta-effect.  psi is the eigenvector.
%
% The y grid is meant to capture the small-scale topographic variation in
% that direction, with no large-scale meridional variation.  The zonal
% variation is large scale relative to the topography, with imposed
% wavenumber k.
%
% see scratch_pde.m for sample input, plotting, and analysis.

% INPUT:
% S: time mean salinity profile (psu)
% T: time mean in situ demperature profile (deg C)
% depth: of model interfaces (m)
% lat: latitude
% k: zonal wavenumber of large scale wave (cpk)
% y: (regular) horizontal grid (m)
% h: topography profile. same size as y. (m)
% optional: number of eigenvalues/eigenvectors to find, and the target
% eigenvalue.  If optional arguments are present, uses eigs (instead of
% eig) to solve the eigenvalue problem.
%
% OUTPUT:
% eigenvecs_mean: 2D array of large-scale eigenvectors.  This is the mean
% of eigenvecs2D in across the horizontal dimension.
% eigenvecs2D: 3D array of eigenvectors.  Will have size 
% [length(depth) length(y) length(w_r)]. The number of eigenvalues is
% variable because I ignore eigenvectors with negative eigenvalues.
% w_r: real part of eigenvals.  (cycles/day)
% w_i: imaginary part of eigencals.  (cycles/day)
% eigenvals: full eigenvalues. (complex)
%
% C. Wortham, 31 Jan. 2020
%
% 23 May 2020: found and fixed error in bottom boundary condition term:
% repeated the wrong element in dh. (Shows up in matrix D.)  Also, was
% making the big dh_matrix wrong.
% 26 May 2020: added option for 2 more inputs: number of eigenvalues to
% find, and the target eigenvalue.  Hoping this will speed things up.
% 29 May 2020: in dealing with negative gradients (static instability), was
% replacing negative gradient with the value below, working my way down.
% But this doesn't work if there are multiple negative values in a row.
% Now work my way up, as I originally intended.
% 17 Feb 2021: Solving the eigenvalue problem with frequency omega as the
% eigenvalue.  Before, had alpha=1/omega as the eigenvalue.  Now also
% output the real and imaginary parts of omega, both sorted in order of
% real part.

%% SORT VARARGIN

optargin = size(varargin,2);

if optargin==0
    solutions = 'full';
elseif optargin==2
    solutions = 'sparse';
    num_vals = varargin{1};
    target_vals = varargin{2};
else
    error('Wrong number of input arguments.  Either pass no optional arguments, or two optional arguments.')
end


%% SETUP

% change to cycles per meter for wavelength
k = k/1000; % (cycles/m)

% global constants (mks units)
global G F RHO0 BETA
G = sw_g(lat,0);
F = sw_f(lat);
RHO0 = 1000;
BETA = 2*7.292e-5*cosd(lat)/6370e3;
% BETA = 0

% other constants (mks units)
Lambda = y(2)-y(1);
Nz = length(depth);
% l is array index in z-direction; l = 1:Nz
Ny = length(y);
% j is array index in y-direction: j = 1:Ny
N = Ny*Nz;
% i is the linear index: i = (j-1)*Nz+l = 1:N

%% BUILD MATRICES

% MATRIX A ~ d^2/dy^2
% diagonal "with fringes" 2nd derivative operator in y
A = (...
    -2*diag(ones(1,N))...
    +diag(ones(1,N-Nz),Nz)...
    +diag(ones(1,N-Nz),-Nz)...
    )/Lambda^2;
% BC's
% TOP       l=1;    i=(j-1)*Nz+1    no change
% BOTTOM    l=Nz;   i=j*Nz          no change
% LEFT      j=1;    i=l
% make a temp matrix with the right form for the first Nz rows
temp = (...
    -2*diag(ones(1,N))...
    +diag(ones(1,N-Nz),Nz)...
    +diag(ones(1,N-(Ny-2)*Nz),(Ny-2)*Nz)...
    )/Lambda^2;
% and replace first Nz rows in A
A(1:Nz,:) = temp(1:Nz,:);
% RIGHT     j=Ny;   i=(Ny-1)*Nz+l
% now do the same with the last Nz rows
temp = (...
    -2*diag(ones(1,N))...
    +diag(ones(1,N-(Ny-2)*Nz),-(Ny-2)*Nz)...
    +diag(ones(1,N-Nz),-Nz)...
    )/Lambda^2;
A(N-Nz+1:end,:) = temp(N-Nz+1:end,:);


% MATRIX B ~ d/dz (F(z)*d/dz)
% tri-diagonal 2nd derivative operator in z
% create some reuseable matrices
n_minus1 = diag(ones(1,N-1),-1);
n = diag(ones(1,N));
n_plus1 = diag(ones(1,N-1),1);
% calculate layer thicknesses
dr = -diff(depth);
dr = ([dr(1) dr] + [dr dr(end)])/2;
% calculate pressure from depth
pres = sw_pres(-depth,lat);
% make sure it's just a vector
dr = squeeze(dr);
if (length(S)~=Nz)
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
% reference pressures
pr_left = ( pres([2 2:Nz]) + pres([1 1:Nz-1]) )/2;
pr_right = ( pres([2:Nz Nz]) + pres([1:Nz-1 Nz-1]) )/2;
% denominator, e.g. rho(n) - rho(n-1), for differential operator 
delta_Rho_left = sw_pden(S([2 2:Nz]),T([2 2:Nz]),pres([2 2:Nz]),pr_left) - sw_pden(S([1 1:Nz-1]),T([1 1:Nz-1]),pres([1 1:Nz-1]),pr_left);
delta_Rho_right = sw_pden(S([2:Nz Nz]),T([2:Nz Nz]),pres([2:Nz Nz]),pr_right) - sw_pden(S([1:Nz-1 Nz-1]),T([1:Nz-1 Nz-1]),pres([1:Nz-1 Nz-1]),pr_right);
% deal with negative gradients
if any(delta_Rho_left<=0)
   if delta_Rho_left(end)<=0
       delta_Rho_left(end) = 1e-5; % make bottom positive
   end
   ind = find(delta_Rho_left<=0);
   for ii=length(ind):-1:1 % work your way up
       delta_Rho_left(ind(ii)) = delta_Rho_left(ind(ii)+1); % work your way up, replacing with next lower value
   end
end
if any(delta_Rho_right<=0)
   if delta_Rho_right(end)<=0
       delta_Rho_right(end) = 1e-5; % make bottom positive
   end
   ind = find(delta_Rho_right<=0);
   for ii=length(ind):-1:1 % work your way up
       delta_Rho_right(ind(ii)) = delta_Rho_right(ind(ii)+1); % work your way up, replacing with next lower value
   end
end

% above is the first step: creating the vector of density differences that
% I need.  Assuming stratification is the same over the whole region (not
% necessary!) I just repeat this vector.
delta_Rho_left_matrix = repmat(delta_Rho_left',[Ny N]);
delta_Rho_right_matrix = repmat(delta_Rho_right',[Ny N]);
dr_matrix = repmat(dr',[Ny N]);
% originally had the following, but realized that the matrix concatenation
% was in the wrong direction:
% delta_Rho_left_matrix = repmat(delta_Rho_left,[N Ny]);
% delta_Rho_right_matrix = repmat(delta_Rho_right,[N Ny]);
% dr_matrix = repmat(dr,[N,Ny]);

% the matrix
B = (F^2*RHO0/G) * dr_matrix.^-1 .* ( (n_minus1 - n)./delta_Rho_left_matrix - (n - n_plus1)./delta_Rho_right_matrix);
% BC's
% I'll do these BC's with a temp matrix again, but here you want to change
% every Nz^th row.
% TOP       l=1;    i=(j-1)*Nz+1
temp = 2*(F^2*RHO0/G) * dr_matrix.^-1 .* (n_plus1 - n)./delta_Rho_right_matrix;
B(1:Nz:end,:) = temp(1:Nz:end,:);
% BOTTOM    l=Nz;   i=j*Nz
temp = 2*(F^2*RHO0/G) * dr_matrix.^-1 .* (n_minus1 - n)./delta_Rho_left_matrix;
B(Nz:Nz:end,:) = temp(Nz:Nz:end,:);
% LEFT      j=1;    i=l             no change
% RIGHT     j=Ny;   i=(Ny-1)*Nz+l   no change


% MATRIX C ~ -4*pi^2*k^2 (from d^2/dx^2)
C = -4*pi^2*k^2 * diag(ones(1,N));
% BC's don't change anything here


% MATRIX D ~ k*BETA (coefficient of eigenvalue)
D = -k*BETA * diag(ones(1,N));
% BC's
% BOTTOM    l=Nz;   i=j*Nz
% the bottom BC includes the effect of topography through the coefficient
% of the eigenvalue omega
% dh = [h(1) diff(h)]; % this was wrong: don't want to repeat the first element of h... fixed below
dh = diff(h);
dh = [dh(1)-dh(end) dh]; % topography is periodic in y, so the extra element in dh should be the difference between the two "ends"
% dh_matrix = repmat(dh',[Nz,N]); % this was wrong way to make the big matrix.
% dh = sqrt(mean(diff(h).^2)) * ones(size(dh)); % just playing around... what happens if I just use the RMS dh/dy?
dh_matrix = zeros(size(D));
for ii=1:Ny
    dh_matrix(ii*Nz,ii*Nz) = dh(ii);
end
topoterm = -(2*F*k/Lambda) .* dr_matrix.^-1 .* dh_matrix .* n;
% now add these elements in D
D(Nz:Nz:end,:) = D(Nz:Nz:end,:) + topoterm(Nz:Nz:end,:);

%% SOLVE EIGENVALUE PROBLEM AND SORT

if strcmp(solutions,'full')
    [eigenvecs,eigenvals] = eig(-D,A+B+C);
else
    [eigenvecs,eigenvals] = eigs(-D,A+B+C,num_vals,target_vals);
end

% % for testing, can ignore A,C,D matrices and make this solve the BTT
% % problem.  Should get same modes as BTTmodes.m.  Also need to change
% % eigenvalue sorting from 'ascend' to 'descend' later on since the
% % interpretation of the eigenvalues changes in this case.  Will get
% % degenerate eigenvalues, with one corresponding to each horizontal grid
% % point.
% [eigenvecs,eigenvals] = eig(B);

% reshape eigenvectors to match physical space.  Just undoing the array to
% linear index switch.  Should preserve mode order.
eigenvecs2D = reshape(eigenvecs,Nz,Ny,[]); % trying to work with eigs

% convert frequency to cycles/day: multiply by 86400 seconds/day
eigenvals = 86400*eigenvals;

% sort output by real part of eigenvalue
eigenvals = diag(eigenvals);
[w_r,ind] = sort(real(eigenvals),'descend');
w_i = imag(eigenvals(ind));
eigenvals = eigenvals(ind);
eigenvecs2D = eigenvecs2D(:,:,ind);

% % % % % comment this out to keep the full solutions.
% % % % % drop any eigenvectors with negative eigenvalues
% % % % ind = find(alpha>=0,1,'first');
% % % % alpha = alpha(ind:end);
% % % % eigenvecs2D = eigenvecs2D(:,:,ind:end);

% redefine N to be number of positive eigenvalues
N = length(w_r);

% % average the eigenvecs2D over horizontal dimension, giving the large scale
% % vertical structure for each mode.
% eigenvecs_mean = squeeze(mean(eigenvecs2D,2));


%% FLIPPING ORIENTATION AND NORMALIZING

% normalize 2D eigenvalues so they integrate to one over the
% domain area
% dr_mat = repmat(dr',[1,Ny,N]);
% AreaInt= -1/depth(end)/y(end)*sum(sum(eigenvecs2D.*conj(eigenvecs2D).*dr_mat,1)*Lambda,2);
% AreaInt_mat = repmat(AreaInt,[Nz,Ny,1]);
% eigenvecs2D = eigenvecs2D./sqrt(AreaInt_mat);
n = 1/(depth(end))/y(end)* trapz(depth,trapz(y,eigenvecs2D.^2,2),1);
n_mat = repmat(n,[size(eigenvecs2D,1),size(eigenvecs2D,2),1]);
eigenvecs2D = eigenvecs2D./sqrt(n_mat);
% set nans (from dividing 0/0) to zero
eigenvecs2D(isnan(eigenvecs2D)) = 0;
% alternative: make maximum to one

% average the eigenvecs2D over horizontal dimension, giving the large scale
% vertical structure for each mode.
eigenvecs_mean = squeeze(mean(eigenvecs2D,2));

% flip mean eigenvecs so that they are positive at the top
% (they are undetermined up to sign)
flips = sign(real(eigenvecs_mean(1,:)));
flips = repmat(flips,[Nz,1]);
eigenvecs_mean = eigenvecs_mean.*flips;

% flip the 2D eigenvecs to match the mean eigenvecs
flips = repmat(flips,[1,1,Ny]);
flips = permute(flips,[1,3,2]);
eigenvecs2D = eigenvecs2D.*flips;


%% dropped code

% took this out... relying on the area normalization above.  Did this as a
% way to ignore extra meridional modes.
% % normalize mean eigenvectors so they integrate to 1 over depth
% dr_mat = repmat(dr',[1,N]);
% VertInt = -1/depth(end)*sum(eigenvecs_mean.*conj(eigenvecs_mean).*dr_mat,1);
% VertInt_mat = repmat(VertInt,[Nz,1]);
% eigenvecs_mean = eigenvecs_mean./sqrt(VertInt_mat);
% % set nans (from dividing 0/0) to zero
% eigenvecs_mean(isnan(eigenvecs_mean)) = 0;
% % alternative: make maximum to one

% % flip eigenvecs so that they alternate sign at the surface (they are
% % undetermined up to sign).  This gets a little convoluted, but I couldn't
% % think of a neater way to do it.
% flips = sign(real(eigenvecs2D(1,1,:)));
% flips = repmat(flips,[Nz,Ny,1]);
% eigenvecs2D = eigenvecs2D.*flips;
% s = (-1).^(2:N+1);
% s = reshape(s,[1 1 N]);
% s = repmat(s,[Nz Ny 1]);
% eigenvecs2D = eigenvecs2D.*s;

% dr_mat = repmat(dr',[1,Ny,N]);
% VertInt= -1/depth(end)*sum(eigenvecs2D.*conj(eigenvecs2D).*dr_mat,1);
% VertInt_mat = repmat(VertInt,[Nz,1,1]);
% eigenvecs2D = eigenvecs2D./sqrt(VertInt_mat);
% % set nans (from dividing 0/0) to zero
% eigenvecs2D(isnan(eigenvecs2D)) = 0;
% % alternative: make maximum to one






