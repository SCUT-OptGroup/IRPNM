clear;

varrho = [0 0.1 0.3 0.6 0.8];

% **************** to fix the random seed **********************
randstate = 100;
randn('state',double(randstate));
rand('state',double(randstate)); 

% get data
loss = 'student';
reg = 'ell1';
n = 262144; 
type = 4;  % type=d/20 and d\in{20,40,60,80}, so you can choose it in {1,2,3,4}, which represents the dynamic range of signal
[xtrue,b,J,k,noise] = getdata(n,type,randstate);
data.b = b;
data.J = J;
m = length(b);
normb = norm(b);
Amap = @(x) pdct(x,1,n,J);
ATmap = @(x) pdct(x,2,n,J);
data.Amap = Amap;
data.ATmap = ATmap; 

% parameter setting
solver = 'DALMQ';
OPTIONS.solver = solver;
OPTIONS.maxiter = 500;
OPTIONS.maxiter_in = 100;
OPTIONS.tol = 1e-8;
OPTIONS.printyes = 1;

model.loss = loss;
model.reg = reg;
nu = 0.25;
model.nu = nu;
t = b./(nu+b.^2);
ATt = pdct(t,2,n,J);
model.lambda = 0.1*2*max(abs(ATt));

pars.eta = 0.9;
pars.b1 = 1.0;

% generate an initial point
x0 = pdct(b,2,n,J); % x0 = ATb

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
        linx1 = 1:scale;
        liny1 = sqrt(sum(dx_mat.^2,1));
        
    elseif pars.varrho == varrho(2)
        xbar = Iter_mat(:,end);
        dx_mat = Iter_mat - xbar;
        scale = size(Iter_mat,2);
        linx2 = 1:scale;
        liny2 = sqrt(sum(dx_mat.^2,1));
        
    elseif pars.varrho == varrho(3)
        xbar = Iter_mat(:,end);
        dx_mat = Iter_mat - xbar;
        scale = size(Iter_mat,2);
        linx3 = 1:scale;
        liny3 = sqrt(sum(dx_mat.^2,1));
        
    elseif pars.varrho == varrho(4)
        xbar = Iter_mat(:,end);
        dx_mat = Iter_mat - xbar;
        scale = size(Iter_mat,2);
        linx4 = 1:scale;
        liny4 = sqrt(sum(dx_mat.^2,1));
        
    elseif pars.varrho == varrho(5)
        xbar = Iter_mat(:,end);
        dx_mat = Iter_mat - xbar;
        scale = size(Iter_mat,2);
        linx5 = 1:scale;
        liny5 = sqrt(sum(dx_mat.^2,1));
    end
end


figure;
semilogy(linx1,liny1,'m-','linewidth',3);
hold on
semilogy(linx2,liny2,'b-','linewidth',3);
semilogy(linx3,liny3,'c-','linewidth',3);
semilogy(linx4,liny4,'r-','linewidth',3);
semilogy(linx5,liny5,'g-','linewidth',3);
xlabel('Number of iterations (with \tau = \rho)');
ylabel('$\log(\|x^k-\overline{x}\|)$','Interpreter','latex');
legend('\rho=0','\rho=0.1','\rho=0.3','\rho=0.6','\rho=0.8','Location','southwest');

