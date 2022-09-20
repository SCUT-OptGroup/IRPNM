%% **************************************************************
%    prox_x = argmin 0.5||z - x||^2 + lambda||w o z||_1
%
%% **************************************************************************

function ML1 = Moreau_L1(x,lambda)

xp = max(0,abs(x)-lambda).*sign(x);

ML1 = 0.5*norm(xp-x)^2+lambda*sum(abs(xp));

end