clear;

varrho = [0   0.1   0.3   0.6   0.8];

% **************** to fix the random seed **********************
randstate = 100;
randn('state',double(randstate));
rand('state',double(randstate)); 

% get data
loss = 'student';
reg = 'ell1';
xtrue = imread('cameraman.tif');
xtrue = double(xtrue);
[p,q] = size(xtrue);
xtrue = xtrue(:);
Amap = @(x) imgaussfilt(x,4,'FilterSize',9,'Padding','symmetric');
Axtrue = Amap(xtrue);
noise = random('T',1,[p*q 1])*1e-3;
b = Axtrue + noise;
Bmap = @(x) phaar(x,1,p,q);
BTmap = @(x) phaar(x,2,p,q);
data.b = b;
data.Amap = @(x) Amap(BTmap(x));
data.ATmap = @(x) Bmap(Amap(x));

% parameter setting
solver = 'DALMQ';
OPTIONS.solver = solver;
OPTIONS.maxiter = 500;
OPTIONS.maxiter_in = 100;
OPTIONS.tol = 3e-8;
OPTIONS.printyes = 1;

model.loss = loss;
model.reg = reg;
model.nu = 1;
model.lambda = 1e-2;

pars.eta = 0.9;
pars.b1 = 1.0;

% generate an initial point
x0 = Bmap(b); 

for i = 1:length(varrho)
    pars.varrho = varrho(i);
    if varrho(i) == 0
        pars.tau = 0
    else
        pars.tau = varrho(i)
    end

    [Iter_mat] = IRPNM_rho(x0, data, model, OPTIONS, pars);
    
    if pars.varrho == varrho(1)
        xbar = Iter_mat(:,end);
        dx_mat = Iter_mat - xbar;
        scale = size(Iter_mat,2);
        linx1_im = 1:scale;
        liny1_im = sqrt(sum(dx_mat.^2,1));
        
    elseif pars.varrho == varrho(2)
        xbar = Iter_mat(:,end);
        dx_mat = Iter_mat - xbar;
        scale = size(Iter_mat,2);
        linx2_im = 1:scale;
        liny2_im = sqrt(sum(dx_mat.^2,1));
        
    elseif pars.varrho == varrho(3)
        xbar = Iter_mat(:,end);
        dx_mat = Iter_mat - xbar;
        scale = size(Iter_mat,2);
        linx3_im = 1:scale;
        liny3_im = sqrt(sum(dx_mat.^2,1));
        
    elseif pars.varrho == varrho(4)
        xbar = Iter_mat(:,end);
        dx_mat = Iter_mat - xbar;
        scale = size(Iter_mat,2);
        linx4_im = 1:scale;
        liny4_im = sqrt(sum(dx_mat.^2,1));
        
    elseif pars.varrho == varrho(5)
        xbar = Iter_mat(:,end);
        dx_mat = Iter_mat - xbar;
        scale = size(Iter_mat,2);
        linx5_im = 1:scale;
        liny5_im = sqrt(sum(dx_mat.^2,1));

    end
end


figure;
semilogy(linx1_im,liny1_im,'m-','linewidth',3);
hold on
semilogy(linx2_im,liny2_im,'b-','linewidth',3);
semilogy(linx3_im,liny3_im,'c-','linewidth',3);
semilogy(linx4_im,liny4_im,'r-','linewidth',3);
semilogy(linx5_im,liny5_im,'g-','linewidth',3);
xlabel('Number of iterations (with \tau = \rho)');
ylabel('$\log(\|x^k-\overline{x}\|)$','Interpreter','latex');
legend('\rho=0','\rho=0.1','\rho=0.3','\rho=0.6','\rho=0.8','Location','southwest');
