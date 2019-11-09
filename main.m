% main.m
% This script when run should compute all values and make all plots
% required by the project.
% To do so, you must fill the functions in the functions/ folder,
% and create scripts in the scripts/ folder to make the required
% plots.

% Add folders to path
addpath('./functions/','./scripts/');

% Add plot defaults
plotDefaults;

%% Exercise 1
T = 1085;
n = 30;
L =1;

[dates,industries,R_m,r_f] = loadStockData1('C:\Users\wc145\Desktop\ECON676\PS4\Problem_Set4.xls');
[Asset_Betas, grsstat,Average_ret] = grs(industries, R_m, r_f);
Rp_mean = mean(industries)';
Rp_std = std(industries)';
Rp_SR = Rp_mean./Rp_std;
R_p = zeros(T,n);
for j = 1:n
    R_p(:,j) = industries(:,j)-r_f;
end
% time series regression
for j = 1:n
    results(j,1) = ols(R_p(:,j),[ones(T,1) R_m]);
end
Cov = zeros(n,n);
Cov_inv = zeros(n,n);
for i = 1:n
    for j = 1:n
        Cov(i,j) = sum(results(i).resid.*results(j).resid) /(T-L-1);
    end
end
Cov_inv = inv(Cov);
mu_m = mean(R_m);
sigma_m = std(R_m);
alpha = zeros(n,1);
for i = 1:n
    alpha(i,1) = results(i).beta(1);
end
W = alpha'*Cov_inv*alpha/(1+mu_m^2/sigma_m^2);
W_normal = T*(T-n-L)*W/(n*(T-L-1));
p = fcdf(W_normal,n,T-n-L,'upper');
%% Exercise 2
T = 1079;
n = 10;
L =1;

[dates,past,R_m,r_f] = loadStockData2('C:\Users\wc145\Desktop\ECON676\PS4\Problem_Set4.xls');
Rp_mean = mean(past)';
Rp_std = std(past)';
Rp_SR = Rp_mean/Rp_std;
R_p = zeros(T,n);
for j = 1:n
    R_p(:,j) = past(:,j)-r_f;
end
% time series regression
for j = 1:n
    results(j,1) = ols(R_p(:,j),[ones(T,1) R_m]);
end
Cov = zeros(n,n);
Cov_inv = zeros(n,n);
for i = 1:n
    for j = 1:n
        Cov(i,j) = sum(results(i).resid.*results(j).resid) /(T-L-1);
    end
end
Cov_inv = inv(Cov);
mu_m = mean(R_m);
sigma_m = std(R_m);
alpha = zeros(n,1);
for i = 1:n
    alpha(i,1) = results(i).beta(1);
end
W = alpha'*Cov_inv*alpha/(1+mu_m^2/sigma_m^2);
W_normal = T*(T-n-L)*W/(n*(T-L-1));
p = fcdf(W_normal,n,T-n-L,'upper');
%% Exercise 3
T = 1085;
n = 25;
L =1;

[dates,beme,R_m,r_f] = loadStockData3('C:\Users\wc145\Desktop\ECON676\PS4\Problem_Set4.xls');
Rp_mean = mean(beme)';
Rp_std = std(beme)';
Rp_SR = Rp_mean/Rp_std;
R_p = zeros(T,n);
for j = 1:n
    R_p(:,j) = beme(:,j)-r_f;
end
% time series regression
for j = 1:n
    results(j,1) = ols(R_p(:,j),[ones(T,1) R_m]);
end
Cov = zeros(n,n);
Cov_inv = zeros(n,n);
for i = 1:n
    for j = 1:n
        Cov(i,j) = sum(results(i).resid.*results(j).resid) /(T-L-1);
    end
end
Cov_inv = inv(Cov);
mu_m = mean(R_m);
sigma_m = std(R_m);
alpha = zeros(n,1);
for i = 1:n
    alpha(i,1) = results(i).beta(1);
end
W = alpha'*Cov_inv*alpha/(1+mu_m^2/sigma_m^2);
W_normal = T*(T-n-L)*W/(n*(T-L-1));
p = fcdf(W_normal,n,T-n-L,'upper');

%% In sample Tangency portfolios
% sample mean and covariance, TP
T = 1085;
n= 25;
 
[dates,beme,R_m,r_f] = loadStockData3('C:\Users\wc145\Desktop\ECON676\PS4\Problem_Set4.xls');
mean_r = zeros(n,1);
for j = 1:n
    mean_r(j,1) = mean(beme(:,j));
end
Rf = mean(r_f);
CovMatrix_r = cov(beme);
[frontier,MVP,MVP_w, TanP, TanP_w]= mvp(mean_r', CovMatrix_r, Rf, 0, 'TitleString', 'PlotName');
r_tp = beme * TanP_w;
R_tp = r_tp-r_f;
% time series regression for TP
R_p = zeros(T,n);
for j = 1:n
    R_p(:,j) = beme(:,j)-r_f;
end
for j = 1:n
    results(j,1) = ols(R_p(:,j),[ones(T,1) R_tp]);
end
Cov = zeros(n,n);
Cov_inv = zeros(n,n);
for i = 1:n
    for j = 1:n
        Cov(i,j) = sum(results(i).resid.*results(j).resid) /(T-L-1);
    end
end
Cov_inv = inv(Cov);
mu_tp = mean(R_tp);
sigma_tp = std(R_tp);
alpha = zeros(n,1);
for i = 1:n
    alpha(i,1) = results(i).beta(1);
end
W = alpha'*Cov_inv*alpha/(1+mu_tp^2/sigma_tp^2);
W_normal = T*(T-n-L)*W/(n*(T-L-1));
p = fcdf(W_normal,n,T-n-L,'upper');
%% Out of Sample Tangency Portfolio
T = 1085;
n = 25;
L = 1;
% Compute values
[dates,beme,R_m,r_f] = loadStockData4('C:\Users\wc145\Desktop\ECON676\PS4\Problem_Set4.xls');
% odd months in even years and even months in odd years
j = 1; k =1;
bemeA = zeros(543,n);
bemeB = zeros(542,n);
r_f_A = zeros(543,1);
r_f_B = zeros(542,1);
for i = 1:T
    if (mod(floor(beme(i,1)/100)+1,2) && mod(beme(i,1),2))||(mod(floor(beme(i,1)/100),2) && mod(beme(i,1)+1,2))
        bemeA(j,:) = beme(i,2:end);
        r_f_A(j,1)= r_f(i,1);
        j = j+1;
    end
end
for i = 1:T
    if (mod(floor(beme(i,1)/100)+1,2) && mod(beme(i,1)+1,2))||(mod(floor(beme(i,1)/100),2) && mod(beme(i,1),2))
        bemeB(k,:) = beme(i,2:end);
        r_f_B(k,1)= r_f(i,1);
        k = k+1;
    end
end   
mean_r_A = zeros(n,1);
for j = 1:n
    mean_r_A(j,1) = mean(bemeA(:,j));
end
Rf_A = mean(r_f_A);
CovMatrix_r_A = cov(bemeA);
[frontier,MVP,MVP_w, TanP, TanP_w_A]= mvp(mean_r_A', CovMatrix_r_A, Rf_A, 0, 'TitleString', 'PlotName');
mean_r_B = zeros(n,1);
for j = 1:n
    mean_r_B(j,1) = mean(bemeB(:,j));
end
Rf_B = nanmean(r_f_B);
CovMatrix_r_B = cov(bemeB);
[frontier,MVP,MVP_w, TanP, TanP_w_B]= mvp(mean_r_B', CovMatrix_r_B, Rf_B, 0, 'TitleString', 'PlotName');
for i = 1:T
    if (mod(floor(beme(i,1)/100)+1,2) && mod(beme(i,1),2))||(mod(floor(beme(i,1)/100),2) && mod(beme(i,1)+1,2))
        r_tp(i,1) = beme(i,2:end)*TanP_w_B;
    end
end
for i = 1:T
    if (mod(floor(beme(i,1)/100)+1,2) && mod(beme(i,1)+1,2))||(mod(floor(beme(i,1)/100),2) && mod(beme(i,1),2))
        r_tp(i,1) = beme(i,2:end)*TanP_w_A;
    end
end        
R_tp = r_tp-r_f;
% time series regression for TP
R_p = zeros(T,n);
for j = 1:n
    R_p(:,j) = beme(:,j+1)-r_f;
end
for j = 1:n
    results(j,:) = ols(R_p(:,j),[ones(T,1) R_tp]);
end
Cov = zeros(n,n);
Cov_inv = zeros(n,n);
for i = 1:n
    for j = 1:n
        Cov(i,j) = sum(results(i).resid.*results(j).resid) /(T-L-1);
    end
end
Cov_inv = inv(Cov);
mu_tp = mean(R_tp);
sigma_tp = std(R_tp);
alpha = zeros(n,1);
for i = 1:n
    alpha(i,1) = results(i).beta(1);
end
W = alpha'*Cov_inv*alpha/(1+mu_tp^2/sigma_tp^2);
W_normal = T*(T-n-L)*W/(n*(T-L-1));
p = fcdf(W_normal,n,T-n-L,'upper');
%% In sample industry tangency
T = 1085;
n = 25;
L =1;
 
[dates,industries,R_m,r_f] = loadStockData1('C:\Users\wc145\Desktop\ECON676\PS4\Problem_Set4.xls');
[dates,beme,R_m,r_f] = loadStockData3('C:\Users\wc145\Desktop\ECON676\PS4\Problem_Set4.xls');
mean_r = zeros(30,1);
for j = 1:30
    mean_r(j,1) = mean(industries(:,j));
end
Rf = mean(r_f);
CovMatrix_r = cov(industries);
[frontier,MVP,MVP_w, TanP, TanP_w]= mvp(mean_r', CovMatrix_r, Rf, 0, 'TitleString', 'PlotName');
r_tp = industries * TanP_w;
R_tp = r_tp-r_f;
% regress on beme
R_p = zeros(T,n);
for j = 1:n
    R_p(:,j) = beme(:,j)-r_f;
end
for j = 1:n
    results(j,1) = ols(R_p(:,j),[ones(T,1) R_tp]);
end
Cov = zeros(n,n);
Cov_inv = zeros(n,n);
for i = 1:n
    for j = 1:n
        Cov(i,j) = sum(results(i).resid.*results(j).resid) /(T-L-1);
    end
end
Cov_inv = inv(Cov);
mu_tp = mean(R_tp);
sigma_tp = std(R_tp);
alpha = zeros(n,1);
for i = 1:n
    alpha(i,1) = results(i).beta(1);
end
W = alpha'*Cov_inv*alpha/(1+mu_tp^2/sigma_tp^2);
W_normal = T*(T-n-L)*W/(n*(T-L-1));
p = fcdf(W_normal,n,T-n-L,'upper');
%% Insample past return tangency
T = 1079;
n = 10;
L =1;

[dates,past,R_m,r_f] = loadStockData2('C:\Users\wc145\Desktop\ECON676\PS4\Problem_Set4.xls');
[dates,beme,R_m,r_f] = loadStockData3('C:\Users\wc145\Desktop\ECON676\PS4\Problem_Set4.xls');
mean_r = zeros(n,1);
for j = 1:n
    mean_r(j,1) = mean(past(:,j));
end
Rf = mean(r_f(7:end,1));
CovMatrix_r = cov(past);
[frontier,MVP,MVP_w, TanP, TanP_w]= mvp(mean_r', CovMatrix_r, Rf, 0, 'TitleString', 'PlotName');
r_tp = past * TanP_w;
R_tp = r_tp-r_f(7:end,1);
% regress on beme
R_p = zeros(T,n);
for j = 1:n
    R_p(:,j) = beme(7:end,j)-r_f(7:end,1);
end
for j = 1:n
    results(j,1) = ols(R_p(:,j),[ones(T,1) R_tp]);
end
Cov = zeros(n,n);
Cov_inv = zeros(n,n);
for i = 1:n
    for j = 1:n
        Cov(i,j) = sum(results(i).resid.*results(j).resid) /(T-L-1);
    end
end
Cov_inv = inv(Cov);
mu_tp = mean(R_tp);
sigma_tp = std(R_tp);
alpha = zeros(n,1);
for i = 1:n
    alpha(i,1) = results(i).beta(1);
end
W = alpha'*Cov_inv*alpha/(1+mu_tp^2/sigma_tp^2);
W_normal = T*(T-n-L)*W/(n*(T-L-1));
p = fcdf(W_normal,n,T-n-L,'upper');
%% STA715 HW2 Q6 Simulation
T = 100;
n = 50;
w = zeros(T,1);
p =  zeros(n*n+1,1);
% k experiments
for k = 1:T
    for i = 1:n
        for j = 1:n
            z(i,j) = normrnd(0,1);
            if z(i,j)<= i*j/(4*n*n)
                w(k,1) = w(k,1) + 1;
            end
        end
    end
end
for i =0:n^2
    for k = 1:T
        if w(k,1) == i
            p(i,1) = p(i,1)+1/(n*n);
        end
    end
end


lambda = 0;
for i = 1:n
    for j = 1:n
        lambda = lambda + i*j/(4*n*n);
    end
end
x = 0:n*n ;     
plot(x,p);
hold on;
plot(x,poisspdf(x,lambda));

%% SML predictibility
T = 1085;
n = 25;
L =1;

[dates,beme,R_m,r_f] = loadStockData3('C:\Users\wc145\Desktop\ECON676\PS4\Problem_Set4.xls');
s = zeros(T,1);
for i = 1:T   
    s(i,1) = sum(beme(i,1:5));   
end
b = zeros(T,1);
for i = 1:T   
    b(i,1) = sum(beme(i,21:25));   
end
r = s./5 - b./5;        
results= ols(r(2:end,1),[ones(T-1,1) r(1:end-1)]);       

%% HML predictibility
T = 1085;
n = 25;
L =1;

[dates,beme,R_m,r_f] = loadStockData3('C:\Users\wc145\Desktop\ECON676\PS4\Problem_Set4.xls');
l = zeros(T,1);
for i = 1:T   
    l(i,1) = beme(i,1)+beme(i,6)+beme(i,11)+beme(i,16)+beme(i,21);   
end
h = zeros(T,1);
for i = 1:T   
    h(i,1) = beme(i,5)+beme(i,10)+beme(i,15)+beme(i,20)+beme(i,25);   
end
r = h./5 - l./5;        
results= ols(r(2:end,1),[ones(T-1,1) r(1:end-1)]);       














