% definition of the Sturm-Liouville BVP
% solving for vertical structure of linear Rossby waves
%
% C. Wortham

function [p,dp,q,r,alp,dom,x,wint,wipr]=BTT_uniform(n)
% The problem -d/dx(p(x)*du/dx)+q(x)*u=lambda*r(x)*u, x in dom=[a,b]
%              alp11*u(a)+alp12*u'(a)=0, alp21*u(b)+alp22*u'(b)=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uniform stratification with rho=1, F=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zmin = -5000;
zmax = 0;
dom= [zmin,zmax] ;% working domain; for example dom=[0,1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% don't touch this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,wint]=pd(n,dom);% Gauss quadrature grid and weights
alpha=(dom(2)-dom(1))/2;
beta=(dom(2)+dom(1))/2;
e1=ones(n,1)/2;e2=e1;e2(2)=1;
X=spdiags([e1,e2],[-1,1],n,n);
X=alpha*X+beta*eye(n);% the matrix X for the domain dom

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set density profile and stratification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=  1;% coefficient of -u''; for example p=X;
dp= 0;% coefficient of -u' (in fact p'(x)); for example dp=1;
q= 0   ;% coefficient of u; for example q=0;
r= 1;% inner product weight (matrix); for example r=2*X*(eye(n)-X^2);
wipr= ones(n,1) ;% inner product weight (vector); for example wipr=2*x.*(1-x.^2);
alp= [0,1;0,1] ;% boundary conditions; for example alp=[0,1;1,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% guide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p = rho(z)*F(z)
% dp = dpdz = drhodz*F + rho*dFdz
% q = 0
% r = rho(z)
% wipr = r in vector form
% alp = matrix of boundary conditions




