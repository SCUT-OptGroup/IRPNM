%% ************************************************************
%   Calculate the projection operator onto B = B1 times B2 times...times Bkappa
%   Bi =\{ z\in R^{|J_i|} | ||zi||<=wi}
%
%% **************************************************************
%% Copyright by Shaohua Pan (2017.10.2)

function projz = Proj_ball(z,pars,gamma,lambda)

ng = pars.ng; gdim = pars.gdim;

if (pars.samedim ==1)   %% each group has the same dimension
    
    Z = reshape(z,gdim,ng);
    
    Zcnorm = dot(Z,Z).^(1/2);      %% which is a row vector
    
    ratio_vec = (1-gamma)+gamma*min(lambda./Zcnorm,1);
    
    Zratio = Z.*ratio_vec;      %Z*diag(ratio_vec);
    
    projz = Zratio(:);
      
else
    
    groups = pars.groups;
    
    projz = z;
   
    m0 = 0;
    
    for i = 1:ng
        
        zi = z(groups==i);
        
        mi = size(zi,1);
              
        normzi = norm(zi);
        
        if (normzi>lambda)
            
            projz(m0+1:m0+mi)=((1-gamma)+gamma*min(lambda/normzi,1))*zi;
            
        end
        
        m0 = m0 + mi;
        
    end
end
