%% ***************************************************************
% To compute the function value of h
%% ***************************************************************

function hval = hval(z,mu,lambda)

[prox,ML1] = Prox_L1(-z/mu,lambda/mu);

hval = 1/(2*mu)*norm(z)^2 -mu*ML1; 

end