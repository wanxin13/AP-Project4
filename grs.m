%% grs.m
%% This function does parts a & b from Weller problem set 4
%% a) Calculate Mean, Stdev & Sharpe Ratio of returns
%% b) Estimate Betas & Intercepts; Do GRS test
%% First factor should be RmRf
%% NOTE: ASSUMES NO MISSING DATA

function [Asset_Betas, grsstat,Average_ret]=grs(Asset_Rets, f, Rf)
% Asset_Category, Names, f_Name)

%% Run time-series regressions to get Beta & intercept of each industry,
%% Also get returns and stdev for each asset when there is data

T = size(Asset_Rets,1);
N = size(Asset_Rets,2);
L = size(f,2);

for i = 1:N
    Y = Asset_Rets(:,i) - Rf;
    X = [f];
    [beta,error,sterrbeta,R2,tstat,param,varbeta]=ordleast(Y,X);
    Asset_Betas(i,:) = beta';
    Asset_Beta_stderrs(i,:) = sterrbeta';
    Asset_Beta_tstats(i,:) = tstat';    
    Resids(:,i) = error;
end

Average_ret = mean(Asset_Rets);
Std_ret = std(Asset_Rets);
Sharpe = sqrt(12)*((Average_ret-mean(Rf))./Std_ret); % Note, the Sharpe ratio is annualized

% %% Part a) Calculate Mean, Stdev & Sharpe Ratio of returns
% disp([Asset_Category '(a) Mean Return, Stdev of return & Sharpe Ratio'])
% Col_Heads = {'Mean';'Stdev';'Sharpe'};
% Results = [Average_ret', Std_ret', Sharpe'];
% make_table(Names, Col_Heads, Results,8,2)
% 
% figure
%   subplot(2,1,1)
%     plot(Std_ret, Average_ret, '*')
%     title([Asset_Category '(a) Stdev vs Returns'])
%     xlabel('Standard deviation of realized returns')
%     ylabel('Average Realized Return')
%   subplot(2,1,2)
%     plot(Average_ret, Sharpe, '*')
%     title([Asset_Category '(a) Average Ret vs Sharpe Ratio'])
%     xlabel('Average Realized Return')
%     ylabel('Sharpe Ratio')
%    
% %% Part b) Estimate Betas & Intercepts; Do GRS test;
% disp([Asset_Category '(b) Intercept, Beta with ' f_Name ' & t-statistics'])
% Col_Heads = {'Intercept';'Beta';'t(int)';'t(Beta)'};
% Results = [Asset_Betas, Asset_Beta_tstats];
% make_table(Names, Col_Heads, Results,9,2)

%% Do GRS test;
% disp([Asset_Category '(b) GRS test for alphas in model with ' f_Name])
T = size(Asset_Rets,1);
N = size(Asset_Rets,2);
L = size(f,2);

mu = mean(f);
sigma=Resids'*Resids/(T-L-1);
omega=cov(f);

alpha=Asset_Betas(:,1);

grsstat=(T/N)*((T-N-L)/(T-L-1))*(alpha'*inv(sigma)*alpha)/(1+mu'*inv(omega)*mu)
grsdof=[N T-N-L]
grsp=1-fcdf(grsstat,grsdof(1),grsdof(2))
grsc5=finv(0.95,grsdof(1),grsdof(2));
grsc1=finv(0.99,grsdof(1),grsdof(2));

% %% Supplemental stuff
% figure
% subplot(2,1,1)
%     plot(sort(Asset_Betas(:,2)),'*')
%     title([Asset_Category ':SUPPLEMENT: Betas (sorted) with ' f_Name]);
% subplot(2,1,2)
%     plot(sort(Asset_Betas(:,1)),'*')
%     title([Asset_Category ':SUPPLEMENT: Intercepts (sorted) in model with ' f_Name]);
