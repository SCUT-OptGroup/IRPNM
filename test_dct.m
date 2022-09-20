clear;

addpath(genpath('DALMQ'));
rmpath(genpath('DALMQ_GP'));

loss = 'student';
reg = 'ell1';

n = 262144; 
tol = 1e-5;
nu = 0.25;

type = 1; % type=d/20 and d\in{20,40,60,80}, so you can choose it in {1,2,3,4}, which represents the dynamic range of signal
c_lam = 0.1; % choose c_lambda in {0.1,0.01}

nrun = 10;

iter_array1 = size(nrun,1); time_array1 = size(nrun,1); 
fval_array1 = size(nrun,1); resi_array1 = size(nrun,1);

for j = 1: nrun
    
    % **************** to fix the random seed **********************
    randstate = j*100
    randn('state',double(randstate));
    rand('state',double(randstate));
    
    % **************** get data **********************
    [xtrue,b,J,k,noise] = getdata(n,type,randstate);
    data.b = b;
    data.J = J;
    m = length(b);
    normb = norm(b);
    Amap = @(x) pdct(x,1,n,J);
    ATmap = @(x) pdct(x,2,n,J);
    data.Amap = Amap;
    data.ATmap = ATmap;

    t = b./(nu+b.^2);
    ATt = pdct(t,2,n,J);
    lambda = c_lam*2*max(abs(ATt));
    
    model.loss = loss;
    model.reg = reg;
    model.nu = nu;
    model.lambda = lambda;
    
    % generate an initial point
    x0 = pdct(b,2,n,J); % x0 = ATb

    % **************** parameter setting of IRPNM **********************
    solver = 'DALMQ';
    OPTIONS.solver = solver;
    OPTIONS.maxiter = 1000;
    OPTIONS.maxiter_in = 100;
    OPTIONS.tol = tol;
    OPTIONS.printyes = 1;
    
    pars.eta = 0.9;
    pars.b1 = 1.0;
    pars.varrho = 0.45;
    pars.tau = 0.45;
    
    [xopt,Fopt,resi,iter,ttime] = IRPNM(x0, data, model, OPTIONS, pars);
    
    iter_array1(j) = iter-1;
    time_array1(j) = ttime;
    fval_array1(j) = Fopt;
    resi_array1(j) = resi;

end

iter1 = mean(iter_array1);  time1 = mean(time_array1);
fval1 = mean(fval_array1);  resi1 = mean(resi_array1);


fprintf('---------------------------------------------------------------------------\n');
fprintf('   Algorithm  |  iter  |  CPU time  |     Obj.    |   Residual  |   Rate    \n');
fprintf('---------------------------------------------------------------------------\n');
fprintf('%12s: \t %3.1f \t  %3.1f \t  %3.4f \t  %3.2e \n', 'IRPNM', iter1, time1, fval1, resi1);
