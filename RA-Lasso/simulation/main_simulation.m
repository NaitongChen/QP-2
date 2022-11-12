cvx_setup

n = 100;  % number of observations
p = 400;  % number of variables
CASE = 5;  % 5 scenarios for error distribution
LAMAX = 1.5;  % lambda_max
AMAX = 1.0;  % alpha_max
LENGTH_la = 30;  % number of lambdas used in validation
LENGTH_a = 15;  % number of alphas used in validation
N = 100;  % number of validation datasets
beta_0 = [3*ones(1,20), zeros(1,380)]';  % true beta

%% Validation for RA-Lasso
XV = csvread('x_ho_v.csv',1);
YV = csvread('y_ho_v.csv',1);
lamb = linspace(0.0001,LAMAX,LENGTH_la);
alp = linspace(0.0001,AMAX,LENGTH_a);

Lambda.Alpha = zeros(CASE,2);  % stores optimal pair of (lambda, alpha) for 6 scenarios

% 1st scenario: N(0,4)
loss1 = zeros(LENGTH_la,LENGTH_a);
for k=1:N
    for i=1:LENGTH_la
        for j=1:LENGTH_a
                [k,i,j]
                [betah,~] = L1PenHuber(YV((n*(k-1)+1):(n*k),1), XV((n*(k-1)+1):(n*k),:), lamb(i), alp(j));
                betah = betah .* (abs(betah) > 1e-04);
                loss1(i,j) = loss1(i,j) + norm(betah-beta_0);
             end
    end
end
[q1,q2,minvalue]=find(loss1==min(loss1(:)));
Lambda.Alpha(1,1) = lamb(q1(1));
Lambda.Alpha(1,2) = alp(q2(1));

% 2nd scenario: 2t_3
loss2 = zeros(LENGTH_la,LENGTH_a);
for k=1:N
    for i=1:LENGTH_la
        for j=1:LENGTH_a
                [k,i,j]
                [betah,~] = L1PenHuber(YV((n*(k-1)+1):(n*k),2), XV((n*(k-1)+1):(n*k),:), lamb(i), alp(j));
                betah = betah .* (abs(betah) > 1e-04);
                loss2(i,j) = loss2(i,j) + norm(betah-beta_0);
             end
    end
end
[q1,q2,minvalue]=find(loss2==min(loss2(:)));
Lambda.Alpha(2,1) = lamb(q1(1));
Lambda.Alpha(2,2) = alp(q2(1));

% 3rd scenario: MixN
loss3 = zeros(LENGTH_la,LENGTH_a);
for k=1:N
    for i=1:LENGTH_la
        for j=1:LENGTH_a
                [k,i,j]
                [betah,~] = L1PenHuber(YV((n*(k-1)+1):(n*k),3), XV((n*(k-1)+1):(n*k),:), lamb(i), alp(j));
                betah = betah .* (abs(betah) > 1e-04);
                loss3(i,j) = loss3(i,j) + norm(betah-beta_0);
             end
    end
end
[q1,q2,minvalue]=find(loss3==min(loss3(:)));
Lambda.Alpha(3,1) = lamb(q1(1));
Lambda.Alpha(3,2) = alp(q2(1));

% 4th scenario: LogNormal
loss4 = zeros(LENGTH_la,LENGTH_a);
for k=1:N
    for i=1:LENGTH_la
        for j=1:LENGTH_a
                [k,i,j]
                [betah,~] = L1PenHuber(YV((n*(k-1)+1):(n*k),4), XV((n*(k-1)+1):(n*k),:), lamb(i), alp(j));
                betah = betah .* (abs(betah) > 1e-04);
                loss4(i,j) = loss4(i,j) + norm(betah-beta_0);
             end
    end
end
[q1,q2,minvalue]=find(loss4==min(loss4(:)));
Lambda.Alpha(4,1) = lamb(q1(1));
Lambda.Alpha(4,2) = alp(q2(1));

% 5th scenario: Weibull
loss5 = zeros(LENGTH_la,LENGTH_a);
for k=1:N
    for i=1:LENGTH_la
        for j=1:LENGTH_a
                [k,i,j]
                [betah,~] = L1PenHuber(YV((n*(k-1)+1):(n*k),5), XV((n*(k-1)+1):(n*k),:), lamb(i), alp(j));
                betah = betah .* (abs(betah) > 1e-04);
                loss5(i,j) = loss5(i,j) + norm(betah-beta_0);
             end
    end
end
[q1,q2,minvalue]=find(loss5==min(loss5(:)));
Lambda.Alpha(5,1) = lamb(q1(1));
Lambda.Alpha(5,2) = alp(q2(1));



%% RA-Lasso
K=100;
l2loss = zeros(K,CASE);
l1loss = zeros(K,CASE);
FP = zeros(K,CASE);
FN = zeros(K,CASE);

X = csvread('x_ho.csv',1);  %% "x_ho.csv" for homoscedastic model, "x_he.csv" for heteroscedastic model
Y = csvread('y_ho.csv',1);

for k=1:K
    Betar = zeros(p,CASE);
    for j = 1:CASE
        [k,j]
        [Betar(:,j),~] = L1PenHuber(Y((n*(k-1)+1):(n*k),j), X((n*(k-1)+1):(n*k),:), Lambda.Alpha(j,1), Lambda.Alpha(j,2));
        Betar(:,j) = Betar(:,j) .* (abs(Betar(:,j)) > 1e-04);
        l2loss(k,j) = norm(Betar(:,j)-beta_0);
        l1loss(k,j) = norm(Betar(:,j)-beta_0,1);
        temp = Betar(:,j);
        FP(k,j) = sum(temp(find(beta_0==0)) ~= 0);
        FN(k,j) = sum(temp(find(beta_0~=0)) == 0);    
    end
end

roundn([mean(l2loss);mean(l1loss);mean(FP);mean(FN)],-2)  %% summarizes the numerical results


