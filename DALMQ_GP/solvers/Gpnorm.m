%% **************************************************************
%  filename: Group_norm
%  
%  To calculate G(x)=(||x_{J_1}||,...,||x_{J_kappa}||)'
%
%% ***************************************************************

function [Gx] = Gpnorm(x,pars)

ng = pars.ng;   gdim = pars.gdim;

if (pars.samedim ==1)    %% each group has the same dimension
    
    X = reshape(x,gdim,ng);
    
    Xcnorm = sum(X.*X,1).^(1/2);
    
    Gx = Xcnorm';
    
else
    
    dimgv = pars.dimgv;
    
    Gx = Gnorm(x,dimgv);
    
end

end

