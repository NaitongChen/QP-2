%Does L-1 Regularization (Robust)
%By A. Emre Barut 2/27/2011

function [beta,MAD] = L1PenL2(Y, X, lambda)

[n,p]=size(X);

cvx_begin
    cvx_quiet(true)
    variable betaHat(p,1)
    minimize sum_square(Y-X*betaHat)+n*lambda*norm(betaHat,1)
cvx_end

beta=betaHat;
MAD=norm(Y-X*betaHat,1)/n;
    
end
