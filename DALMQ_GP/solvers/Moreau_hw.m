%% **************************************************************
%  filename: Moreau_hw
%
%% ****************** the envelope of L21 **************************
%
%  to compute the Moreau-envelop of hk(x) = ||w o x||_21
%  
%  ehk = 0.5*gamma*norm(xp-x)^2 + ||w o xp||_21
%
%  xp = argmin{0.5*gamma*||y-x||^2 + ||w o y||_21}
%
%% ****************************************************************

function ehk = Moreau_hw(z,gamma,pars,lambda)

ng = pars.ng;  gdim = pars.gdim;

lam_gam = lambda/gamma;

if (pars.samedim ==1)   %% each group has the same dimension
    
    Z = reshape(z,gdim,ng);
    
    Zcnorm = dot(Z,Z).^(1/2);
   
    ratio_vec = max(1-lam_gam./Zcnorm,0);
    
    Zpmat = Z.*ratio_vec;
    
   % zp = Zpmat(:);
     
    ehk = gamma*(0.5*norm(Zpmat-Z,'fro')^2+lam_gam*sum(Zcnorm));  
    
else
    n = size(z,1);
    
    dimgv = pars.dimgv;
    
    zgnorm = Gnorm(z,dimgv);
    
    zp = zeros(n,1);
    
    pidx = find(zgnorm>lambdaw/gam);
    
    for i=1:length(pidx)
        
        kk = pidx(i);
        
        if kk==1
            
            jj = 0;
        else
            jj = sum(dimgv(1:kk-1));
        end
        
        Ji = [jj+1:jj+dimgv(kk)];

        tempzi = (1-lam_gam/zgnorm(kk))*z(Ji);
        
        zp(Ji) = tempzi;  
    end
      
     ehk = 0.5*gam*norm(zp-z)^2+lam_gam*sum(zgnorm);  
 
end    

