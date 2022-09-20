%% ***************************************************************
%  To compute the function and gradient values of Phi_j  
%
%% ***************************************************************

function [fval,proxh_u,zeta,xnz_ind] = fgrad_DALMQ(xi,u,lambda,sigma,mu)

%gamma = mu*sigma/(1+mu*sigma);

proxh_u = Proxh(u,1/(1+mu*sigma),lambda);

zeta = u - proxh_u;

xnz_ind = zeta~=0;  

hval = 1/(2*mu)*norm(proxh_u)^2 -mu*Moreau_L1(-proxh_u/mu,lambda/mu);

fval = 0.5*norm(xi)^2 + (sigma/2)*norm(zeta)^2 + hval; %(proxh_u,mu,lambda);

end