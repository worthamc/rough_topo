function Lz = BTTmatrix(S,T,pres,dr)
% VERTICAL_STRUCTURE_FUNCTION The vertical component of the
%  quasi-geostrophic PV operator.
%  (NOTE: Still needs to be multiplied by f^2 rho0 / g)
%
% Lz = BTTmatrix(S,T,pres,dr)
%
%   Lz: a finite-difference matrix operator
%   S: salinity (psu)
%   T: in situ temperature (deg C)
%   pres: pressure
%   dr: spacing between centers of grid points
%
%   The stratification is assumed to be stable at all levels.
%
% Modified from R. Abernathey's vertical_structure_function.m
% C. Wortham, 24 Nov. 2018
%
% 29 May 2020: in dealing with negative gradients (static instability), was
% replacing negative gradient with the value below, working my way down.
% But this doesn't work if there are multiple negative values in a row.
% Now work my way up, as I originally intended.


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

Lz = zeros(N);

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
	disp('negative density gradients found')
   if delta_Rho_left(end)<=0
       delta_Rho_left(end) = 1e-8; % make bottom positive
   end
   ind = find(delta_Rho_left<=0);
   for ii=length(ind):-1:1 % work your way up
       delta_Rho_left(ind(ii)) = delta_Rho_left(ind(ii)+1); % work your way up, replacing with next lower value
   end
end
if any(delta_Rho_right<=0)
	disp('negative density gradients found')
   if delta_Rho_right(end)<=0
       delta_Rho_right(end) = 1e-8; % make bottom positive
   end
   ind = find(delta_Rho_right<=0);
   for ii=length(ind):-1:1 % work your way up
       delta_Rho_right(ind(ii)) = delta_Rho_right(ind(ii)+1); % work your way up, replacing with next lower value
   end
end

% make matrices
delta_Rho_left_matrix = repmat(delta_Rho_left',[1 N]);
delta_Rho_right_matrix = repmat(delta_Rho_right',[1 N]);

% delta_Rho_left = Rho - Rho([1 1:N-1]);
% delta_Rho_right = Rho([2:N N]) - Rho;
% delta_Rho_left_matrix = repmat(delta_Rho_left',[1 N]);
% delta_Rho_right_matrix = repmat(delta_Rho_right',[1 N]);

dr_matrix = repmat(dr',[1 N]);

Lz = dr_matrix.^-1 .* ( (n_minus1 - n)./delta_Rho_left_matrix - (n - n_plus1)./delta_Rho_right_matrix);

% boundary conditions
Lz([1 end],:) = zeros(2,N);
Lz(1,1:2) = 2*[-1 1] / (dr(1) * delta_Rho_right(1));
Lz(N,[N-1 N]) = 2*[1 -1] / (dr(N) * delta_Rho_left(N));
% % alternative: phi(z=-H) = 0
% Lz(N,[N-1 N]) = [0 0];

end

