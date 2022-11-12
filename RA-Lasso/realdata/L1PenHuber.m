%Does L-1 Regularization (Robust) using Huber loss
%By Yuyan 05/16/2014

function [beta,MAD] = L1PenHuber(Y, X, lambda, alpha)

[n,p]=size(X);

cvx_begin
    cvx_quiet(true)
    variable betaHat(p,1)
    minimize sum(huber(Y-X*betaHat,1/alpha))+n*lambda*norm(betaHat,1)     %% M = 1\alpha in huber func
cvx_end

beta=betaHat;
MAD=norm(Y-X*betaHat,1)/n;
    
end