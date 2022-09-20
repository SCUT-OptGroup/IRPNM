%% **************************************************************
%    prox_x = argmin 0.5||z - x||^2 + lambda||w o z||_1
%
%% **************************************************************************

function [xp,ML1] = Prox_L1(x,lambdaw)

xp = max(0,abs(x)-lambdaw).*sign(x);

if nargout>=2
    ML1 = 0.5*norm(xp-x)^2+lambdaw*sum(abs(xp));
end
end

