function Iter_mat = IRPNM_rho(x0, data, model, OPTIONS, pars)

% This function uses inexact regularized proximal Newton to solve
% l1-regularized regression problem.

%% get parameters

if isfield(OPTIONS,'maxiter');      maxiter     = OPTIONS.maxiter;       end
if isfield(OPTIONS,'maxiter_in');   maxiter_in  = OPTIONS.maxiter_in;    end
if isfield(OPTIONS,'tol');          tol         = OPTIONS.tol;           end
if isfield(OPTIONS,'printyes');     printyes    = OPTIONS.printyes;      end

n = size(x0,1);

loss = model.loss;

switch loss
    case 'student'
        lambda = model.lambda;
        b = data.b;
        m = size(b,1);
        dataAmap = data.Amap;
        dataATmap = data.ATmap;
end

b1 = pars.b1;
varrho = pars.varrho;
tau = pars.tau;
eta = pars.eta;
info.rhoflag = 1;

if varrho ==0
    tau = 0;
    info.rhoflag = 0;
end

if (printyes)
    fprintf('\n *****************************************************');
    fprintf('******************************************');
    fprintf('\n \t IRPNM for solving %s regularization %s regression with %s innersolver',model.reg,model.loss,OPTIONS.solver);
    fprintf('\n *****************************************************');
    fprintf('*******************************************');
    fprintf('\n  iter   iter_in    time     residual        dnorm         resi_in         Fval           Lambda         sigma');
end

%% Initialization

delta = 1e-4;

beta = 0.1;

%% parameters for DALMQ

OPTIONS_DALMQ.tol = tol;

OPTIONS_DALMQ.maxiter = maxiter_in;

OPTIONS_DALMQ.printyes = 0;

%% main loop

iter = 1;

total_in = 0;

xi = zeros(m,1);

Iter_mat = [];

% Compute function value, gradient and residual at x0

x = x0; 

tstart = clock;

Iter_mat =[Iter_mat  x];

Ax = dataAmap(x);

[Fval,gval,grad,Dx] = Fun_grad(x,Ax,b,data,model,pars,lambda);

resi = norm(x-Prox(x-grad,model));

resi_array(iter) = resi;

b2 = min(1e-2/max(1,resi),delta);    %% this is robust!

while (iter<=maxiter)
    
    if resi<=tol
        
        fprintf('\n IRPNM achieves the required accuracy with residual %3.2e and total inner iteration %3.0d\n',resi,total_in);
        
        Iter_mat =[Iter_mat  x];
        
        return;
    end
    
    switch loss

        case 'student'
            % Compute Lambda and mu at x where the hessian of f takes the form of A'*diag(Dx)*A
            
            Lambda = -b1*min(0,min(Dx));
            
             if Lambda>1e-4    %% corresponds to non-strong stationary point
                
               info.flag_SN = 0;
               
               OPTIONS_DALMQ.sigma = 10;
               
            else
                
               info.flag_SN = 1;
               
               OPTIONS_DALMQ.sigma = 200;   %[100,500]             
            end
            
            mu = b2*resi^varrho;
            
            % **************** solve the subproblem with DALMQ *************
            % initialize Ainput = D_{+}^(1/2)A and binput
            
            Dxp = Dx + Lambda;
            
            adax = dataATmap(Dxp.*Ax);
            
            binput = adax + mu*x - grad;
            
            normb = norm(binput);
            
            sqrtDxp = Dxp.^(1/2);
            
            Amap = @(x) sqrtDxp.*dataAmap(x);
            
            ATmap = @(x) dataATmap(sqrtDxp.*x);
            
            OPTIONS_DALMQ.normb = normb;
            
%             OPTIONS_DALMQ.epsin = max(tol,eta*min(resi,resi^(1+tau)));
            OPTIONS_DALMQ.epsin = eta*min(resi,resi^(1+tau));
            
            info.qks = 0.5*(dot(x,adax)+mu*norm(x)^2)-dot(grad,x);
            
            info.gval = gval;
            
            [y,xi,iter_in,resi_in,sigma] = DALMQ(xi,ATmap(xi),-x,Amap,ATmap,binput,OPTIONS_DALMQ,pars,info,lambda,mu);

    end
    
    total_in = total_in + iter_in;
    %===============================================================================================
    
    Ay = dataAmap(y);
    
    Fval_y = Fun_grad(y,Ay,b,data,model,pars,lambda);
    
    d = y - x;  Ad = Ay - Ax;
    
    dnorm = norm(d);
    
    if (dnorm <= tol)
        
        fprintf('\n IRPNM achieves the required accuracy with dnorm %3.2e and total inner iteration %3.0d\n',dnorm,total_in);
        
        xopt = y;  Fopt = Fval_y;
        
        Iter_mat =[Iter_mat  x];

        return;
    end
    
    %% ******************* Perform line search ********************
    
    ls = 0;
    
    descent = delta*mu*dnorm^2;
    
    xnew = y;  Axnew = Ay;  Fval_new = Fval_y;
    
    while (Fval_new > Fval - beta^(ls)*descent)
        
        ls = ls + 1;
        
        xnew = x + beta^(ls)*d;
        
        Axnew = Ax + beta^(ls)*Ad;
        
        Fval_new = Fun_grad(xnew,Axnew,b,data,model,pars,lambda);
    end
    
    % Update iterate
    
    if (Fval_y < Fval_new)
        x = y; Ax = Ay;
    else
        x = xnew; Ax = Axnew;
    end
    
    ttime = etime(clock,tstart);
    
    Iter_mat =[Iter_mat  x];

    [Fval,gval,grad,Dx] = Fun_grad(x,Ax,b,data,model,pars,lambda);
    
    resi = norm(x-Prox(x-grad,model));
    
    if (printyes)

        fprintf('\n %3.0d     %3.0d       %3.2f      %3.2e       %3.2e      %3.2e      %3.4f      %3.3e     %3.2e',iter,iter_in,ttime,resi,dnorm,resi_in,Fval,Lambda,sigma)
        
    end
    
    iter = iter + 1;
    
    resi_array(iter) = resi;
    
end

Iter_mat =[Iter_mat  x];

end




    
    
        
        
    
        
    
    
