%% ***************************************************************
% The conjugate_gradient method for the system of linear equations
% 
%        Ax = b    with Ax = x + gamma*AmapJ(V*AtmapJ(x))
%
% x0£º the starting point
% 
% res£ºthe residual of Ax = b at x0, i.e., res = b- A(x0) 
%% ***************************************************************
function [x,Atx,iter,solve_ok] = cg_DALMQ(Amap,ATmap,uvec,b,pars,lambda,gamma_mu,res,tol,maxit,nr)
 
m = size(b,1); resnrm =[];

if ~exist('tol','var'); tol=1e-2*norm(b); end
 
if ~exist('maxit','var'); maxit = 10; end
 
solve_ok = 1;
 
tiny = 1.0e-16;
 
stagnate_check = 20;

d = zeros(m,1);   

Ad = zeros(m,1);  Atd = zeros(nr,1);

%% *************** Initialization part **************************

x = zeros(m,1);     % such a starting point is crucial !!!

Atx = Atd;

r = res;    z = r;  
 
err = norm(r);

resnrm(1) = err;  minres = err;

tau_old = err;  rho_old = err^2;
 
theta_old = 0;

%% ******************* Main Loop ********************************

for iter = 1:maxit
    
    Atz = ATmap(z); 

    [Az,sigma] = Jacobian_DALMQ(Amap,uvec,z,Atz,pars,lambda,gamma_mu);
 
    if (abs(sigma)<tiny)   %% in this case z=0 since A is positive definite
        
        solve_ok = -1;
        
        return;        
    else
        
        alfa = rho_old/sigma;
        
        r = r - alfa*Az;    

    end
    
    norm_r = norm(r);
    
    theta = norm_r/tau_old;
    
    c = 1/sqrt(1+theta^2);  
    
    tau = norm_r*c;     
    
    gam = (c*theta_old)^2;
    
    eta = alfa*c^2;  
   
    d =  gam*d + eta*z;
    
    x = x + d;
    
 %%--------------------- stopping conditions  ---------------------------
    
    Ad = gam*Ad + eta*Az;
    
    Atd = gam*Atd + eta*Atz;
    
    Atx = Atx + Atd;
    
    res = res - Ad;
    
    err = norm(res);
    
    resnrm(iter+1) = err;
    
    [~,indictor] = find(err < minres);
    
    minres(indictor) = err(indictor);
    
    if (err < tol)
       
        return; 
    end
    
    if (iter > stagnate_check) && (iter > 10)
        
        ratio = resnrm(iter-9:iter+1)./resnrm(iter-10:iter);
        
        if (min(ratio)>0.997) && (max(ratio)<1.003)
            
            solve_ok = -1;
            
            return;
        end
    end
 %%-----------------------------------------------------------------------
    
    rho = norm_r^2;
    
    beta = rho/rho_old;
    
    z = r + beta*z;
            
    rho_old = rho;
    
    tau_old = tau;
    
    theta_old = theta; 
    
end
 
if (iter == maxit); solve_ok = -2; end
 
end    

%%%%%%%%%%%%%%%%%%%%%% End of conjugate_gradient.m  %%%%%%%%%%%%%%%%%%%%%%%   