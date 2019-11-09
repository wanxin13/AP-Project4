function results=ols(y,x)
% PURPOSE: least-squares regression 
%---------------------------------------------------
% USAGE: results = ols(y,x)
% where: y = dependent variable vector    (nobs x 1)
%        x = independent variables matrix (nobs x nvar)
%---------------------------------------------------
% RETURNS: a structure
%        results.meth  = 'ols'
%        results.beta  = bhat     (nvar x 1)
%        results.tstat = t-stats  (nvar x 1)
%        results.bstd  = std deviations for bhat (nvar x 1)
%        results.yhat  = yhat     (nobs x 1)
%        results.resid = residuals (nobs x 1)
%        results.sige  = e'*e/(n-k)   scalar
%        results.rsqr  = rsquared     scalar
%        results.rbar  = rbar-squared scalar
%        results.dw    = Durbin-Watson Statistic
%        results.nobs  = nobs
%        results.nvar  = nvars
%        results.y     = y data vector (nobs x 1)
%        results.bint  = (nvar x2 ) vector with 95% confidence intervals on beta
%---------------------------------------------------
% SEE ALSO: prt(results), plt(results)
%---------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com
%
% Barry Dillon (CICG Equity)
% added the 95% confidence intervals on bhat

if (nargin ~= 2); error('Wrong # of arguments to ols'); 
else
 [nobs nvar] = size(x); [nobs2 junk] = size(y);
 if (nobs ~= nobs2); error('x and y must have same # obs in ols'); 
 end;
end;

results.meth = 'ols';
results.y = y;
results.nobs = nobs;
results.nvar = nvar;

if nobs < 10000
  [q r] = qr(x,0);
  xpxi = (r'*r)\eye(nvar);
else % use Cholesky for very large problems
  xpxi = (x'*x)\eye(nvar);
end;

results.beta = xpxi*(x'*y);
results.yhat = x*results.beta;
results.resid = y - results.yhat;
sigu = results.resid'*results.resid;
results.sige = sigu/(nobs-nvar);
tmp = (results.sige)*(diag(xpxi));
sigb=sqrt(tmp);
results.bstd = sigb;
tcrit=-tdis_inv(.025,nobs);
results.bint=[results.beta-tcrit.*sigb, results.beta+tcrit.*sigb];
results.tstat = results.beta./(sqrt(tmp));
ym = y - mean(y);
rsqr1 = sigu;
rsqr2 = ym'*ym;
results.rsqr = 1.0 - rsqr1/rsqr2; % r-squared
rsqr1 = rsqr1/(nobs-nvar);
rsqr2 = rsqr2/(nobs-1.0);
if rsqr2 ~= 0
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
else
    results.rbar = results.rsqr;
end;
ediff = results.resid(2:nobs) - results.resid(1:nobs-1);
results.dw = (ediff'*ediff)/sigu; % durbin-watson


function x = beta_inv(p, a, b)
% PURPOSE: inverse of the cdf (quantile) of the beta(a,b) distribution
%--------------------------------------------------------------
% USAGE: x = beta_inv(p,a,b)
% where:   p = vector of probabilities
%          a = beta distribution parameter, a = scalar
%          b = beta distribution parameter  b = scalar
% NOTE: mean [beta(a,b)] = a/(a+b), variance = ab/((a+b)*(a+b)*(a+b+1))
%--------------------------------------------------------------
% RETURNS: x at each element of p for the beta(a,b) distribution
%--------------------------------------------------------------
% SEE ALSO: beta_d, beta_pdf, beta_inv, beta_rnd
%--------------------------------------------------------------

%       Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg
% documentation modified by LeSage to
% match the format of the econometrics toolbox

if (nargin ~= 3)
    error('Wrong # of arguments to beta_inv');
end
 
if any(any((a<=0)|(b<=0)))
   error('beta_inv parameter a or b is nonpositive');
end
if any(any(abs(2*p-1)>1))
   error('beta_inv: A probability should be 0<=p<=1');
end

x = a ./ (a+b);
dx = 1;
while any(any(abs(dx)>256*eps*max(x,1)))
   dx = (betainc(x,a,b) - p) ./ beta_pdf(x,a,b);
   x = x - dx;
   x = x + (dx - x) / 2 .* (x<0);
end

function pdf = beta_pdf(x, a, b)
% PURPOSE: pdf of the beta(a,b) distribution
%--------------------------------------------------------------
% USAGE: pdf = beta_pdf(x,a,b)
% where:   x = vector of components
%          a = beta distribution parameter, a = scalar
%          b = beta distribution parameter  b = scalar
% NOTE: mean[(beta(a,b)] = a/(a+b), variance = ab/((a+b)*(a+b)*(a+b+1))
%--------------------------------------------------------------
% RETURNS: pdf at each element of x of the beta(a,b) distribution
%--------------------------------------------------------------
% SEE ALSO: beta_d, beta_pdf, beta_inv, beta_rnd
%--------------------------------------------------------------

%       Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg
% documentation modified by LeSage to
% match the format of the econometrics toolbox
  

if (nargin ~=3)
    error('Wrong # of arguments to beta_pdf');
end

if any(any((a<=0)|(b<=0)))
   error('Parameter a or b is nonpositive');
end

I = find((x<0)|(x>1));

pdf = x.^(a-1) .* (1-x).^(b-1) ./ beta(a,b);
pdf(I) = 0*I;


function x = tdis_inv (p, a)
% PURPOSE: returns the inverse (quantile) at x of the t(n) distribution
%---------------------------------------------------
% USAGE: x = tdis_inv(p,n)
% where: p = a vector of probabilities 
%        n = a scalar dof parameter
%---------------------------------------------------
% RETURNS:
%        a vector of tinv at each element of x of the t(n) distribution      
% --------------------------------------------------
% SEE ALSO: tdis_cdf, tdis_rnd, tdis_pdf, tdis_prb
%---------------------------------------------------

%       Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg
   

  if (nargin ~= 2)
    error ('Wrong # of arguments to tdis_inv');
  end

s = p<0.5; 
p = p + (1-2*p).*s;
p = 1-(2*(1-p));
x = beta_inv(p,1/2,a/2);
x = x.*a./((1-x));
x = (1-2*s).*sqrt(x);

