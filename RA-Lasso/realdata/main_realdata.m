%% Real Data Analysis
X = csvread('x_real.csv',1);
Y = csvread('y_real.csv',1);
n = size(X,1);
p = size(X,2);
X = [ones(n,1) X];  % add 1's as the 1st column of X corresponding to the intercept term

%% Five-fold cross-validation
fold_index = zeros(5,24);
fold_index(1,:) = linspace(1,24,24);
fold_index(2,:) = linspace(25,48,24);
fold_index(3,:) = linspace(49,72,24);
fold_index(4,:) = linspace(73,96,24);
fold_index(5,:) = linspace(96,119,24);  

nlambda = 10; % number of lambdas used in validation
nalpha = 10; % number of alphas used in validation
lambda1 = linspace(0.1,5,nlambda);
lambda2 = linspace(0.1,5,nlambda);
lambda3 = linspace(0.1,5,nlambda);
rss1 = zeros(1,nlambda);  % RSS for Lasso
rss2 = zeros(1,nlambda);  % RSS for R-Lasso
rss3 = zeros(nalpha,nlambda);  % RSS for RA-Lasso

%% Cross-validation for Lasso
nlambda = 100;
lambda1 = exp(linspace(log(0.001),log(50),nlambda));
rss1 = zeros(1,nlambda);
for j = 1:nlambda
    for k = 1:5
       [j,k]
       x = removerows(X,'ind',fold_index(k,:));
       y = removerows(Y,'ind',fold_index(k,:));
       [beta1,~] = L1PenL2(y, x, lambda1(j));
       beta1 = beta1 .* (abs(beta1) > 1e-03);
       rss1(j) = rss1(j) + norm(Y(fold_index(k,:))-X(fold_index(k,:),:)*beta1,1)/length(fold_index(1,:));
    end
end
rss1 = rss1/5;
rss1
[~,minind] = min(rss1);
lambda1(minind)   
plot(log(lambda1),rss1)   %% CV-plot
xlabel('log(lambda)','FontSize',16);
ylabel('L1-loss','FontSize',16);
lambda_lasso = 9.7052;  %% pick the optimal lambda such that the CV-curve changes most dramatically
line([log(lambda_lasso),log(lambda_lasso)],[0,rss1(find(abs(lambda1-lambda_lasso)<1e-4))],'LineStyle','--');


%% Cross-validation for R-Lasso
nlambda = 100;
lambda2 = exp(linspace(log(0.05),log(0.5),nlambda));
rss2 = zeros(1,nlambda);
for j = 1:nlambda
    for k = 1:5
       [j,k]
       x = removerows(X,'ind',fold_index(k,:));
       y = removerows(Y,'ind',fold_index(k,:));
       [beta2,~] = L1PenL1Single(y, x, lambda2(j));
       beta2 = beta2 .* (abs(beta2) > 1e-04);
       rss2(j) = rss2(j) + norm(Y(fold_index(k,:))-X(fold_index(k,:),:)*beta2,1)/length(fold_index(1,:));
    end
end
rss2 = rss2/5;
rss2
[~,minind] = min(rss2);
lambda2(minind)    %% 0.0207
plot(log(lambda2(40:185)),rss2(40:185))
xlabel('log(lambda)','FontSize',16);
ylabel('L1-loss','FontSize',16);
lambda_rlasso = 0.1972;
line([log(lambda_rlasso),log(lambda_rlasso)],[0.4,rss2(find(abs(lambda2-lambda_rlasso)<1e-4))],'LineStyle','--');


%% Cross-validation for RA-Lasso
nlambda = 100;
nalpha = 5;
lambda3 = linspace(0.01,0.05,nlambda);
alpha = linspace(4,6,nalpha);
rss3 = zeros(nlambda,nalpha);
for nl = 1:nlambda
    for na = 1:1
        [nl,na]
        for k = 1:5
           x = removerows(X,'ind',fold_index(k,:));
           y = removerows(Y,'ind',fold_index(k,:));
           beta3 = zeros(p, 1);
           [beta3,~] = L1PenHuber(y, x, lambda3(nl),12.5750);
           beta3 = beta3 .* (abs(beta3) > 1e-04);
           rss3(nl,na) = rss3(nl,na) + norm(Y(fold_index(k,:))-X(fold_index(k,:),:)*beta3,1)/length(fold_index(1,:));
        end
    end
end       
rss3 = rss3/5;
rss3
[minl,mina] = find(rss3==min(rss3(:)))
lambda3(minl)  
alpha(mina)
plot(log(lambda3(15:100)),rss3(15:100,1))
xlabel('log(lambda)','FontSize',16);
ylabel('L1-loss','FontSize',16);
lambda_ra = 0.0225;
alpha_ra = 5;
line([log(lambda_ra),log(lambda_ra)],[30,rss3(find(abs(lambda3-lambda_ra)<1e-4),1)],'LineStyle','--');


%% Compute Lasso, R-Lasso and RA-Lasso estimates
[Beta1,~] = L1PenL2(Y, X, lambda_lasso);
 Beta1 = Beta1 .* (abs(Beta1) > 1e-03);
 [Beta2,~] = L1PenL1Single(Y, X, lambda_rlasso);  
 Beta2 = Beta2 .* (abs(Beta2) > 1e-03);
 [Beta3,~] = L1PenHuber(Y, X, lambda_ra, alpha_ra); 
 Beta3 = Beta3 .* (abs(Beta3) > 1e-03);
 
 find(Beta1~=0)
 find(Beta2~=0)
 find(Beta3~=0)
 length(find(Beta1~=0))
 length(find(Beta2~=0))
 length(find(Beta3~=0))  %% find which and how many variables are selected by the 3 methods
 
 error1 = Y-X*Beta1;
 qqplot(error1(find(Y-X*Beta1 < 2.5)))
 title('')
 error2 = Y-X*Beta2;
 qqplot(error2(find(Y-X*Beta1 < 2.5)))
  title('')
 error3 = Y-X*Beta3;
 qqplot(error3(find(Y-X*Beta1 < 2.5)))
  title('')
 
lambda_lasso = 9.7052;
lambda_rlasso = 0.0854;
lambda_ra = 0.0128;
alpha_ra = 5;

%% Use Cross-validation to compute the errors for Lasso, R-Lasso and RA-Lasso
N = 100;
L2_lasso = zeros(N,1);
L1_lasso = zeros(N,1); 
L2_rlasso = zeros(N,1);
L1_rlasso = zeros(N,1);
L2_ra = zeros(N,1);
L1_ra = zeros(N,1);
for k = 1:N
    k
    r=randperm(119);
    test_index = r(1:24);
    all_index = 1:119;
    train_index = all_index(~ismember(all_index,test_index));
    train_x = X(train_index,:);
    train_y = Y(train_index);
    test_x = X(test_index,:);
    test_y = Y(test_index);
    [beta1,~] = L1PenL2(train_y, train_x, lambda_lasso);  
    beta1 = beta1 .* (abs(beta1) > 1e-03);
    [beta2,~] = L1PenL1Single(train_y, train_x, lambda_rlasso);  
    beta2 = beta2 .* (abs(beta2) > 1e-03);
    [beta3,~] = L1PenHuber(train_y, train_x, lambda_ra, alpha_ra);
    beta3 = beta3 .* (abs(beta3) > 1e-03);
    L2_lasso(k) = norm(test_y - test_x*beta1,2)^2/length(test_index);
    L1_lasso(k) = norm(test_y - test_x*beta1,1)/length(test_index);
    L2_rlasso(k) = norm(test_y - test_x*beta2,2)^2/length(test_index);
    L1_rlasso(k) = norm(test_y - test_x*beta2,1)/length(test_index);
    L2_ra(k) = norm(test_y - test_x*beta3,2)^2/length(test_index);
    L1_ra(k) = norm(test_y - test_x*beta3,1)/length(test_index);
end
L2_lasso
L1_lasso
L2_rlasso
L1_rlasso
L2_ra
L1_ra

     
[median(L2_lasso), median(L2_rlasso), median(L2_ra)]   
[median(L1_lasso), median(L1_rlasso), median(L1_ra)]   


roundn([median(L2_lasso), median(L2_rlasso), median(L2_ra)],-2)  
roundn([median(L1_lasso), median(L1_rlasso), median(L1_ra)],-2)  

roundn([mean(L2_lasso), mean(L2_rlasso), mean(L2_ra)],-2) 
roundn([mean(L1_lasso), mean(L1_rlasso), mean(L1_ra)],-2)  

lasso	= [rss1; lambda1];
rlasso	= [rss2; lambda2];
ralasso = [rss3 lambda3']';

csvwrite('lasso.csv', lasso);
csvwrite('rlasso.csv', rlasso);
csvwrite('ralasso.csv', ralasso);

L2loss = [L2_lasso, L2_rlasso, L2_ra];
L1loss = [L1_lasso, L1_rlasso, L1_ra];
csvwrite('L2loss.csv',L2loss);
csvwrite('L1loss.csv',L1loss);

%%CC cleaner



 
    
    
    
    


