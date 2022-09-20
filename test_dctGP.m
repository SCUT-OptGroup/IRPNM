clear;

addpath(genpath('DALMQ_GP'));
rmpath(genpath('DALMQ'));

loss = 'student';

reg = 'L21';

n = 262144; 

ng = 1024;                  % the number of groups

gdim = floor(n/ng);       % dimension of every group

groups = ceil([1:n]'/gdim);    % the structure of groups

%% ***************** to formulate AmapJ and AtmapJ ***************

gnact = 16;   %[16  64  128]   % active groups

group_flag = 1;     % 1 for sequential group structure, i.e. groups = ceil([1:p]'/gd);

pars.ng = ng;

pars.gdim = floor(n/ng);

pars.sn = floor(n/8); 

pars.groups = ceil([1:n]'/gdim);    % the structure of groups  

pars.samedim = 1;        % 1 for the same dimension; otherwise 0

if pars.samedim==0
    
    pars.dimgv = gdim*ones(ng,1);   %[16; gdim*ones(ng-1,1)];
    
    n = sum(pars.dimgv);
end  

tol = 1e-5;

nu = 0.2;

type = 3;  % type=d/20 and d\in{60,80}, so you can choose it in {3,4}, which represents the dynamic range of signal

nrun = 10;

iter_array1 = size(nrun,1); time_array1 = size(nrun,1); 
fval_array1 = size(nrun,1); resi_array1 = size(nrun,1);

for j = 1:nrun
    
   % ********************** to fix the random seed **********************
    randstate = j*100
    randn('state',double(randstate));
    rand('state',double(randstate));
    
  % **************************** get data ********************************
    [xtrue,b,J,k,noise] = getdataGP(n,type,nu,gnact,pars,randstate);
    
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
    lambda = 0.1*2*max(abs(ATt));
    
    model.loss = loss;
    model.reg = reg;
    model.nu = nu;
    model.lambda = lambda;
    
    % generate an initial point
    x0 = pdct(b,2,n,J);   %x0 = ATb

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
    [xopt,Fopt,resi_array,iter,ttime] = IRPNM_GP(x0, data, model, OPTIONS, pars);
    
    iter_array1(j) = iter-1;
    time_array1(j) = ttime;
    fval_array1(j) = Fopt;
    resi_array1(j) = resi_array(end);

end

iter1 = mean(iter_array1);  time1 = mean(time_array1);
fval1 = mean(fval_array1);  resi1 = mean(resi_array1);


fprintf('---------------------------------------------------------------------------\n');
fprintf('   Algorithm  |  iter  |    Obj.   |    Residual   |   CPU time  \n');
fprintf('---------------------------------------------------------------------------\n');
fprintf('%12s: \t %3.1f \t  %3.4f \t  %3.2e \t  %3.2f  \n', 'IRPNM', iter1, fval1, resi1,time1);
