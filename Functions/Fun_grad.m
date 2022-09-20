%% ***************************************************************
% This function returns the objective value and gradient 
%% ***************************************************************
function [Fval,gval,grad,dd] = Fun_grad(x,Ax,b,data,model,pars,lambda)

loss = model.loss;

switch loss
 
    case 'student'
        nu = model.nu;
        Axb = Ax-b;
        Axb_sq = Axb.^2;
        fval = sum(log(1+Axb_sq/nu));
        if nargout>=2
            grad = 2*data.ATmap(Axb./(nu+Axb_sq));
            dd = 2*(nu-Axb_sq)./(nu+Axb_sq).^2;
        end

end

reg = model.reg;

switch reg
    
    case 'ell1'
        gval = lambda*norm(x,1);
    case 'L21'
        gval = lambda*sum(Gpnorm(x,pars));

end

Fval = fval+gval;
end
