%% *************************************************************************
%  filename: Prox_c21norm 
%
%% ************************************************************************
%% This code is to solve the following simple convex optimization 
%
%  Y = argmin {0.5*||Y-G||^2 + gamma*sum(||Yi||_2)}
%  
%  where Yi is the ith column of the matrix Y.  
%
%    y = argmin{0.5*||y-g||^2+ gamma*||y||_2}, 
%
%    y* = max(1-gamma/||g||,0)g  if ||g|| is not equal to 0
%
%% *************************************************************************
%% Copyright by Ruyu Li and Shaohua Pan, 2002/08/23

function proxz = Prox_c21(z,pars,gamma)

if norm(z)== 0
    
    proxz = 0;
    
    return;
else
    ng = pars.ng; gdim = pars.gdim;
    
    if (pars.samedim ==1)   %% each group has the same dimension
        
        proxz_mat = zeros(gdim,ng);
        
        Z = reshape(z,gdim,ng);
        
        Zcnorm = sum(Z.*Z,1).^(1/2);      %% which is a row vector
        
        proxz_mat(:,Zcnorm>0) = Z(:,Zcnorm>0).*max(1-gamma./Zcnorm(Zcnorm>0),0);
        
        proxz = proxz_mat(:);
        
    end
end   