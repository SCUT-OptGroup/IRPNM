%% ***************************************************************
%  Consider the L_1-norm regularized problem
%
%  min  0.5*||z||^2 - <b,x> + lambda*||x||_1 + 0.5*mu||x||^2: 
%  s.t. Ax-z=0,                                                  (*) 
%
%  The dual of the problem in (*) takes the following formulation 
%
%  min 0.5*||xi||^2 + h(zeta)
%  s.t. A'*xi-zeta-b=0,                                  (**)
%  where h(z) = (1/2mu)||z||^2 - e_{g/mu}(-z/mu)
%
% This code is the ALM+SNewton for solving the above dual problem (**)

%% ***************************************************************
%% Copyright by Ruyu Liu and Shaohua Pan (2022.4.14)

function [xvar,xi,iter,resi_in,sigma] = DALMQ(xi,Atxi,x,Amap,ATmap,b,OPTIONS,pars,info,lambda,mu)

if isfield(OPTIONS,'epsin');       epsin     = OPTIONS.epsin;    end
if isfield(OPTIONS,'printyes');    printyes  = OPTIONS.printyes;   end
if isfield(OPTIONS,'maxiter');     maxiter   = OPTIONS.maxiter;    end
if isfield(OPTIONS,'normb');       normb     = OPTIONS.normb;      end
if isfield(OPTIONS,'sigma');       sigma     = OPTIONS.sigma;      end

sigma_max = 1e+5; 

sigma_min = 1e-5;

rfac = 1.2; 

total_SN = 0;

prim_win = 0;  dual_win = 0;  

if (printyes)
    fprintf('\n *****************************************************');
    fprintf('******************************************');
    fprintf('\n ************** ALM for the constrained weighted L1-norm regularized problem ****************');
    fprintf('\n ****************************************************');
    fprintf('*******************************************');
    fprintf('\n  iter itSN       resi_in          distq        sigma   time');
end

%% ********** parameters for the semismooth newton method ************

OPTIONS_SNDir.printyes = 0;

if info.flag_SN==1
    
    OPTIONS_SNDir.maxiter = 10;
else
    OPTIONS_SNDir.maxiter = 3;
end

OPTIONS_SNDir.tol = epsin;

OPTIONS_SNDir.normb = normb;
   
%% ************** Initialization part for ALMSN ******************

tstart = clock;

%% ******************************* Main Loop ******************************

for iter = 1:maxiter
    
    %********************** solve the ALM subproblem ******************
    
    [xi,Atxi,proju,x,Ax,itSN] = SNCG_DALMQ(xi,Atxi,Amap,ATmap,b,OPTIONS_SNDir,x,pars,lambda,sigma,mu);
    
    total_SN = total_SN + itSN;
    
 %***********************  update the Lagrange multiplier *******************

    xvar = -x;  As = -Ax; 
    
    pobj = 0.5*norm(As)^2-b'*xvar+lambda*sum(Gpnorm(xvar,pars))+(mu/2)*norm(xvar)^2;
    
    AtAs = ATmap(As);
    
% *************************** stopping condition ***************************
    
    if info.rhoflag == 1
        
        Gxvar = AtAs+mu*xvar;
        
        prox_res = xvar - Prox_c21(xvar+b-Gxvar,pars,lambda);
        
        resi_in = norm(prox_res);
    else
        resi_in = norm(AtAs-proju-b);
    end
    
    primfeas = norm(xi-As)/(1+normb);
    
    dualfeas = norm(Atxi-proju-b)/(1+norm(proju));
    
    if (primfeas<dualfeas)
        
        prim_win = prim_win+1;
    else
        dual_win = dual_win+1;
    end
    
    ttime = etime(clock,tstart);
    
    if (printyes) 
        
        fprintf('\n %3d    %3d       %3.5e       %3.2f  %3.2f',iter,itSN,resi_in,sigma,ttime);
    end

    if (resi_in<=epsin)&&(pobj+info.qks<=info.gval)
        
        return;
    end
    
 %% *************** update the parameters and variables **************
    
    sigma_update_iter = 1;
    
    if mod(iter,sigma_update_iter) ==0
        if (prim_win>max(1,dual_win))
            prim_win = 0.0;
            sigma = min(sigma_max,sigma*rfac);
            
        elseif (dual_win>max(1,prim_win))
            dual_win = 0.0;
            sigma = max(sigma_min,sigma/rfac);
        end
   end 

end
 
%     if iter <=10
%        sigma_update_iter = 1;
%       elseif iter <=100
%            sigma_update_iter = 3;           
%       elseif iter <=200
%            sigma_update_iter = 5;
%       else
%            sigma_update_iter = 10;
%       end