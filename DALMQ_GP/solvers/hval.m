%% ***************************************************************
% To compute the function value of h
%% ***************************************************************

function hval = hval(z,mu,pars,lambda)

[~,ML1] = Moreau_hw(-z/mu,mu,pars,lambda);

hval = 1/(2*mu)*norm(z)^2 - ML1; 

end