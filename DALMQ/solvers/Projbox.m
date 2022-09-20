%% ************************************************************
%   Calculate the projection operator onto the set 
%   B =\{ z\in R^p} | |z_i|<=omega_i*gamma }
%% **************************************************************

function projz = Projbox(z,lambda)

projz = z;

% zabs = abs(z);
% 
% tempw = lambda*sign(z(zabs>lambda));
% 
% projz(zabs>lambda) = tempw; 

projz(abs(z)>lambda) = lambda*sign(z(abs(z)>lambda));

