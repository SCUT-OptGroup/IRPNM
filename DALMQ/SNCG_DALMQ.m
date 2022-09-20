%% ******** Semismooth Newton Method for solving the system ********
%
%    nabla Phi_k(xi)=0
%
%% ***************************************************************
%% Copyright by Shaohua Pan and Ruyu Liu (2022.4.14)

function [xi,Atxi,proju,negx,Anegx,iter]= SNCG_DALMQ(xi,Atxi,Amap,ATmap,b,OPTIONS,x,lambda,sigma,mu)

%% ******************* Initialization of parameter *********************

if isfield(OPTIONS,'tol');             tolSN      = OPTIONS.tol;         end
if isfield(OPTIONS,'printyes');        printyes   = OPTIONS.printyes;    end
if isfield(OPTIONS,'maxiter');         maxiter    = OPTIONS.maxiter;     end
if isfield(OPTIONS,'normb');           normb      = OPTIONS.normb;       end

sigma_mu = sigma/(1+mu*sigma);

cg_tol = min(1e-2,100*tolSN);

cg_max = 5;

if  (printyes)
    fprintf('\n *****************************************************');
    fprintf('******************************************');
    fprintf('\n \t   Semismooth Newton method for solving the problem of min phik(xi)');
    fprintf('\n ****************************************************');
    fprintf('*******************************************');
    fprintf('\n  iter cg_iter   lstep     grad_res         obj     time');
end

%% *************** the parameters for line search ***************

lsopts.c1 = 1e-4;

lsopts.c2 = 0.9;

lsopts.stepop = 1;

lsopts.printyes = 0;

nr = size(x,1);

%% ******************** Main Loop *******************************

xsigma = x/sigma;

u = Atxi + xsigma - b;

[fval,proju,zeta,xnz_ind] = fgrad_DALMQ(xi,u,lambda,sigma,mu);

gradPsi = xi + sigma*Amap(zeta);

if abs(fval)>1.0e+5
    ls_tol = 1e-4;    
elseif abs(fval)>1.0e+2    
    ls_tol = 1e-5;
else
    ls_tol = 1e-6;
end

res = -gradPsi;

norm_res = norm(res);

for iter = 1:maxiter
    if (norm_res<=tolSN) 
        
        negx = sigma*zeta;
        
        Anegx = gradPsi - xi;
        
        return;
    end
        
    cg_tol = min(cg_tol,norm_res^(1.2));
    
    %% ******************* to calculate Newton direction **********************
    
    [dir,Atdir,cg_iter] = cg_DALMQ(Amap,ATmap,xnz_ind,res,sigma_mu,res,cg_tol,cg_max,nr);
        
    [xi,Atxi,proju,zeta,xnz_ind,fval,gradPsi,lstep]=...
        Lswolfe_DALMQ(xi,Atxi,fval,gradPsi,dir,Atdir,Amap,b,xsigma,lambda,sigma,mu,lsopts,ls_tol);
    
    res = -gradPsi;
    
    norm_res = norm(res);
    
    if (printyes)        
        fprintf('\n %3.0d      %3.0d      %3.2e     %3.2e    %5.4e   %3.2e  %3.2e',iter,cg_iter,lstep,norm_res,fval,sigma,cg_tol);
    end  
   
            
end

negx = sigma*zeta;

Anegx = gradPsi - xi;


