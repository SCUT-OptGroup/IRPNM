function [xtrue,b,J,k,noise,A] = getdata(n,i,randstate)

randn('state',double(randstate));
rand('state',double(randstate)); 

 m = floor(n/8);
 
 % Generate xtrue
 k = floor(n/40);
 xtrue = zeros(n,1);
 ind = randperm(n,k);
 eta1 = randsrc(k,1);
 eta2 = rand(k,1);
 xtrue(ind) = eta1.*10.^(i*eta2);
 
 % Generate b = A*xtrue + noise, where A*xtrue = dct(xtrue)
 J = randperm(n,m);
 temp = dct(xtrue);
 Axtrue = temp(J);
 noise = random('T',4,[m 1]);
 b = Axtrue + 1e-1*noise;
% b = Axtrue;
 
 % calculate the dct matrix
 if n <= 28500
     B = dctmtx(n);
     A = B(J,:);
 end
     
 