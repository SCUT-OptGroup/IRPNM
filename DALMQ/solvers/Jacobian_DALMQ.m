%% ***************************************************************
%  Given vector d, calculate 
%
%  Ad = d + gamma/(1+mu*gamma)*Anz(Anzt(d))
%  
%  sigma = <d,Ad> = ||d||^2 + gamma/(1+mu*gamma)*||Anzt(d)||^2
%
%% **************************************************************

function [Ad,sigma] = Jacobian_DALMQ(Amap,Anz_ind,d,Atd,gamma_mu)

Anztd = Atd.*Anz_ind;

Ad = d + gamma_mu*Amap(Anztd);

sigma = norm(d)^2 + gamma_mu*norm(Anztd)^2;

end

