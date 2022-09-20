function px = Prox(x,model,pars)
reg = model.reg;
switch reg
    case 'ell1'
        lambda = model.lambda;
        px = Prox_L1(x,lambda);
    case 'simplex'
        px = Proj_simplex(x,1);
    case 'Haar'
        p = model.p;  
        q = model.q;
        lambda = model.lambda;
        Bx = phaar(x,1,p,q);
        prox_Bx = Prox_L1(Bx,lambda);
        px = phaar(prox_Bx,2,p,q);
    case 'L21'
        lambda = model.lambda;
        px = Prox_c21(x,pars,lambda);
end
end