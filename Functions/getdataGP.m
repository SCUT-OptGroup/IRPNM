function [xtrue,b,J,k,noise,A] = getdataGP(n,i,nu,gnact,pars,randstate)

randn('state',double(randstate));
rand('state',double(randstate)); 

 m = pars.sn;
 
 gdim = pars.gdim;
 
 % Generate xtrue
 k = gnact*gdim;
 xtrue = zeros(n,1);
 ind = randperm(n,k);
 eta1 = randsrc(k,1);
 eta2 = rand(k,1);
 xtrue(ind) = eta1.*10.^(i*eta2);
 
 % Generate b = A*xtrue + noise, where A*xtrue = dct(xtrue)
 J = randperm(n,m);
 temp = dct(xtrue);
 Axtrue = temp(J);
 noise = random('T',floor(1/nu),[m 1]);
 b = Axtrue + 1e-1*noise;
% b = Axtrue;
 
 % calculate the dct matrix
 if n <= 28500
     B = dctmtx(n);
     A = B(J,:);
 end
     
 