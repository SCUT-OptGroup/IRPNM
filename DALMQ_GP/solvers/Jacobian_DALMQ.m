%% ***************************************************************

%  Calculate the Clarke generalized Jacobian for the projection Pi_Bi

%  V = I - Jac(Pi_Bi)  and  Vz = V*z;
%% **************************************************************
%%  Copyright by Ruyu Liu and Shaohua Pan (2022/08/22)

function [Az,sigma] = Jacobian_DALMQ(Amap,u,z,Atz,pars,lambda,gamma_mu)

ng = pars.ng;   gdim = pars.gdim;

if (pars.samedim ==1)    %% each group has the same dimension
    
    VZ = zeros(gdim,ng);
    
    U = reshape(u,gdim,ng);
    
    Ucnorm = dot(U,U).^(1/2);    % the column sum of U
    
    Uratio = U.*(1./Ucnorm);
    
    Z = reshape(Atz,gdim,ng);
    
    ZU = dot(Z,Uratio);
    
    alpha = lambda./Ucnorm;   %% alpha is a row vector
    
    alpha_ZU = alpha.*ZU;     %% alpha_ZU is a row vector
    
    if sum(Ucnorm>lambda)>0
        
        tempW1 = Z(:,Ucnorm>lambda).*(1-alpha(Ucnorm>lambda));
        
        tempW2 = Uratio(:,Ucnorm>lambda).*alpha_ZU(Ucnorm>lambda);
        
        VZ(:,Ucnorm>lambda) = tempW1+tempW2;
    end
    
    Az = z + gamma_mu*Amap(VZ(:));
    
    sigma = z'*Az;
    
else
    
    dimgv = pars.dimgv;
    
    ugnorm = Gnorm(u,dimgv);
    
    Jidx = find(ugnorm>=lambdaw);
    
    Vz = zeros(size(Atz,1),1);
    
    for i = 1:length(Jidx)
        
        kk = Jidx(i);
        
        if kk==1
            
            jj = 0;
        else
            jj = sum(dimgv(1:kk-1));
        end
        
        Ji = [jj+1:jj+dimgv(kk)];
        
        ui = u(Ji);  uinorm = norm(ui);
       
        uibar = ui/uinorm;
       
        alpha0 =lambda/uinorm;
        
        Atzi = Atz(Ji);
        
        Vz(Ji)=(1-alpha0)*Atzi + alpha0*(uibar'*Atzi)*uibar;
    end
    
    Az = z + gamma_mu*Amap(Vz);
    
    sigma = z'*Az;
end





