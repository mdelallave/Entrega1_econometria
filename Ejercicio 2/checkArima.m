function check_ar = checkArima(y,pp,D,qq,T)
%% checkArima: function to check ARIMA parameters - ARIMA(p,d,q)
%
%% SYNTAX:
%  check_ar = checkArima(y,pp,qq)
%
%% INPUT:
%  y : Data
%  pp: Maximum for p
%  D : Number of the diferences
%  qq: Maximum for q
%  T : Number of observations in the data set 
%% OUTPUT:
%  ar : matrix of the Bayesian Information Criterion (BIC) of the estimated
%       models where the rows correspond to the AR degree (p) and the
%       columns correspond to the MA degree (q). The smallest value is 
%       the best model.
%
%% EXAMPLE:
%  pp = 4; D = 1; qq = 4; y = data;
%  check_ar = checkArima(y,pp,d,qq)
%

LOGL = zeros(pp,qq); %Initialize
PQ = zeros(pp,qq);
for p = 1:pp
    for q = 1:qq
        model = arima(p,D,q)
        [EstMdl,~,logL] = estimate(model,y,'Display','off');
        LOGL(p,q) = logL;
        PQ(p,q) = p + q;
     end
end
LOGL = reshape(LOGL,pp*qq,1);
PQ = reshape(PQ,pp*qq,1);
[~,bic] = aicbic(LOGL,PQ+1,T);
check_ar = reshape(bic,pp,qq);
