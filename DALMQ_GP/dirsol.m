%% ***************************************************************
%  filename: dirsol
%
%% ****************************************************************
%  Calculate the system: Udir + A*V*A'dir = Res by using Woodbury formula
%
%  U: the Clarke generalized Jacobian of Prox_{1/gamma2}f(zj+u/gamma2)
%
%  V: the Clarke generalized Jacobian of A Prox_{1/gamma1}hk(xj-A'*u/gamma1) 
%
%% **************************************************************

function [dir] = dirsol(u,A,Res,pars,lambdaw,sigma)

[n,p] = size(A);  

gdim = pars.gdim;

%% ********************** Jacobian of Proxhk ***********************

if (pars.samedim ==1)   %% each group has the same dimension
    
    ugnorm = Group_norm(u,pars);
    
    nzidx = find(ugnorm>=lambdaw);
    
    sumAVA = 0;
           
    for i = 1:length(nzidx)
        
        kk = nzidx(i);
        
        Ji = [(kk-1)*gdim+1:kk*gdim];
           
        Ai = A(:,Ji);
                       
        ui = u(Ji);  uinorm = norm(kk);
            
        uibar = ui/uinorm;
            
        Auibar = Ai*uibar;
        
        aki = lambdaw(kk)/uinorm;
                      
        sumAVA = sumAVA + (sigma*(1-aki)*Ai)*Ai'+ (sigma*aki*Auibar)*Auibar';
        
    end
    
    sumAVA = (sumAVA+sumAVA')/2;
    
    dir =(eye(n)+sumAVA)\Res;

else   
    
    dimgv = pars.dimgv;
    
    ugnorm = Gnorm(zeta,dimgv);
    
    nzidx = find(ugnorm>=lambdaw);
    
    sumAVA = zeros(n);
           
    for i = 1:length(nzidx)
        
        kk = nzidx(i);
        
        if kk==1
            
           jj = 0;           
        else   
           jj = sum(dimgv(1:kk-1));           
        end
        
        Ji = [jj+1:jj+dimgv(kk)];
        
        ui = u(Ji);  uinorm = ugnorm(kk); 
        
        uibar = ui/uinorm;
        
        aki = lambdaw(kk)/uinorm;
               
        Ai = A(:,Ji);   
                   
        Auibar = Ai*uibar;
            
        sumAVA = sumAVA + (sigma*(1-aki)*Ai)*Ai'+(aki*Auibar)*Auibar';
        
    end
      
    sumAVA = (sumAVA+sumAVA')/2;
    
    dir =(eye(n)+sumAVA)\Res;
end        
end

