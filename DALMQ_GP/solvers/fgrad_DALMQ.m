%% ***************************************************************
%  To compute the function and gradient values of Phi_j  
%
%% ***************************************************************

function [fval,proxh_u,zeta] = fgrad_DALMQ(xi,u,pars,lambda,sigma,mu)

proxh_u = Proj_ball(u,pars,1/(1+mu*sigma),lambda);

zeta = u - proxh_u;

hval = 1/(2*mu)*norm(proxh_u)^2 - Moreau_hw(-proxh_u/mu,mu,pars,lambda); 

fval = 0.5*norm(xi)^2 +(sigma/2)*norm(zeta)^2 + hval;   %(proxh_u,mu,pars,lambda);

end