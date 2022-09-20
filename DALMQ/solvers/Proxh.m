%% ***************************************************************
%
% Compute the proximal operator of h
%
%% ***************************************************************

function proxv = Proxh(v,gamma,lambda)

tempv = v;

tempv(abs(v)>lambda) = lambda*sign(v(abs(v)>lambda));

proxv = (1-gamma)*v + gamma*tempv;

end
