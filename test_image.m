clear;

addpath(genpath('DALMQ'));
rmpath(genpath('DALMQ_GP'));

loss = 'student';
reg = 'ell1';

maxiter = 1000;
tol = 1e-4;

nrun = 10;

iter_array1 = size(nrun,1); time_array1 = size(nrun,1); 
fval_array1 = size(nrun,1); resi_array1 = size(nrun,1); 
ISNR_array1 = size(nrun,1);

xtrue = imread('cameraman.tif');
xtrue = double(xtrue);
[p,q] = size(xtrue);
xtrue = xtrue(:);
Amap = @(x) imgaussfilt(x,4,'FilterSize',9,'Padding','symmetric');
Axtrue = Amap(xtrue);

for j = 1:nrun
    % **************** to fix the random seed **********************
    randstate = j*100
    randn('state',double(randstate));
    rand('state',double(randstate));
    
    % get data
    noise = random('T',1,[p*q 1])*1e-3;
    b = Axtrue + noise;
    Bmap = @(x) phaar(x,1,p,q);
    BTmap = @(x) phaar(x,2,p,q);
    data.b = b;
    data.Amap = @(x) Amap(BTmap(x));
    data.ATmap = @(x) Bmap(Amap(x));
    
    nu = 1;
    lambda = 1e-2;
    model.loss = loss;
    model.reg = reg;
    model.nu = nu;
    model.p = p;
    model.q = q;
    model.lambda = lambda;

    pars.eta = 0.9;
    pars.b1 = 1.0;
    pars.varrho = 0.45;
    pars.tau = 0.45

    % parameter setting
    solver = 'DALMQ';
    OPTIONS.solver = solver;
    OPTIONS.maxiter = maxiter;
    OPTIONS.maxiter_in = 100;
    OPTIONS.tol = tol;
    OPTIONS.printyes = 1;
    
    x0 = Bmap(b);
    
    [xopt,Fopt,resi,iter,ttime] = IRPNM(x0, data, model,OPTIONS, pars);
    
    xopt = BTmap(xopt);
    obs_IRPNM = uint8(xopt);
    obs_IRPNM = reshape(obs_IRPNM,[p q]);
    subplot(1,3,3),imshow(obs_IRPNM);title('IRPNM');
    
    iter_array1(j) = iter - 1;
    time_array1(j) = ttime;
    fval_array1(j) = Fopt;
    resi_array1(j) = resi;
    ISNR_array1(j)= 10*log10(norm(xtrue-b)^2/norm(xopt-xtrue)^2);
 
end
%% reconstruct the picture and comparision
noisepic = uint8(b);
noisepic = reshape(noisepic,[p q]);
ori = reshape(uint8(xtrue),[p q]);
subplot(1,3,1),imshow(noisepic); title('noisy blurred image');
subplot(1,3,2),imshow(ori);title('original image');

%% display the average result
iter1 = mean(iter_array1);
time1 = mean(time_array1);
fval1 = mean(fval_array1);
resi1 = mean(resi_array1);
ISNR1 = mean(ISNR_array1);


fprintf('---------------------------------------------------------------------------\n');
fprintf('   Algorithm  |  iter  |   Obj.  |   Residual  | CPU time |   ISNR   \n');
fprintf('---------------------------------------------------------------------------\n');
fprintf('%12s: \t %3.1f \t  %3.4f \t  %3.2e \t  %3.2f \t %3.4f \n', 'IRPNM', iter1, fval1, resi1,time1,ISNR1);




% % save image
% figure;imshow(obs_IRPNM,'Border','tight'); saveas(gcf,'./Figures/fig_imIRPNM.eps','epsc');saveas(gcf,'./Figures/fig_imIRPNM.fig')
