clear all
%% Muestra
M = 1.0e4;
T = [20,50,100];
mu = 0;
sigma = 1;
signif_level = [0.01, 0.05, 0.1];
X1 = normrnd(mu,sigma,M,T(1));
X2 = normrnd(mu,sigma,M,T(2));
X3 = normrnd(mu,sigma,M,T(3));

JB_X1= zeros(M,length(signif_level));
JB_X2= zeros(M,length(signif_level));
JB_X3= zeros(M,length(signif_level));

for j = 1:length(signif_level)
    for i = 1:M
        JB_X1(i,j) = jbtest(X1(i,:),signif_level(j));
        JB_X2(i,j) = jbtest(X2(i,:),signif_level(j));
        JB_X3(i,j) = jbtest(X3(i,:),signif_level(j));
    end
end

n_rejects_X1 = sum(JB_X1 == 1,1);
per_rejects_X1 = n_rejects_X1/M*100;
n_rejects_X2 = sum(JB_X2 == 1,1);
per_rejects_X2 = n_rejects_X2/M*100;
n_rejects_X3 = sum(JB_X3 == 1,1);
per_rejects_X3 = n_rejects_X3/M*100;


%%%%% Compute our own JB statistic %%%
% sk_X1 = skewness(X1(1,:),1); % No corregir el sesgo
% kur_X1 = kurtosis(X1(1,:),1); % No corregir el sesgo
% 
% myJB_X1 = T(1)/6*(sk_X1.^2 + 0.25*(kur_X1-3).^2) % Our statistic
% [h,p,jbstat,critval] = jbtest(X1(1,:),0.05) % Matlab's function
% % Results:
% % myJB_X1 =
% % 
% %     1.4387
% % jbstat =
% % 
% %     1.4387
%%%%%

%% Contraste de potencia

Z1 = trnd(3,M,T(1));
Z2 = trnd(3,M,T(2));
Z3 = trnd(3,M,T(3));

JB_Z1= zeros(M,length(signif_level));
JB_Z2= zeros(M,length(signif_level));
JB_Z3= zeros(M,length(signif_level));

for j = 1:length(signif_level)
    for i = 1:M
        JB_Z1(i,j) = jbtest(Z1(i,:),signif_level(j));
        JB_Z2(i,j) = jbtest(Z2(i,:),signif_level(j));
        JB_Z3(i,j) = jbtest(Z3(i,:),signif_level(j));
    end
end

n_rejects_Z1 = sum(JB_Z1 == 1,1);
per_rejects_Z1 = n_rejects_Z1/M*100;
n_rejects_Z2 = sum(JB_Z2 == 1,1);
per_rejects_Z2 = n_rejects_Z2/M*100;
n_rejects_Z3 = sum(JB_Z3 == 1,1);
per_rejects_Z3 = n_rejects_Z3/M*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRUEBAS ULTIMO APARTADO %%%%%%%%%%%%%%%%%%%%%
sample_mean = mean(Z1(1,:)); % provisional, extender a todas las muestras

% Varianzas poblacionales calculadas como la media de las varianzas muestrales
var_pobl_Z1 = mean(var(Z1,0,2));

sk_Z1 = skewness(Z1,0,2);
kur_Z1 = kurtosis(Z1,0,2);
myJB_Z1 = T(1)/6*(sk_Z1.^2 + 0.25*(kur_Z1-3).^2);
critical_v = chi2inv(1-signif_level,2);
test_result = myJB_Z1 < critical_v; 
n_rejected_Z1 = M - sum(test_result);
per_rejected = n_rejected_Z1/M;

% dc para el contraste unilateral. H0:mu=5; H1: mu>5
dc_uni=norminv(1-signif_level(1),0,1);
% dc para el contraste bilateral. H0:mu=5; H1: mu<>5
dc_bi=norminv(1-signif_level(1)/2,0,1);
%Región Crítica para el contraste unilateral
RC_uni=mu+(sigma/(T(1)^0.5))*dc_uni;
%Región Crítica para el contraste bilateral
RC_bi=[mu-(sigma/(T(1)^0.5))*dc_bi,mu+(sigma/(T(1)^0.5))*dc_bi];
%p-value para el contraste unilateral
pvalue_uni=1-normcdf(W0,0,1);
%p-value para el contraste bilateral
pvalue_bi=2*(1-normcdf(W0,0,1));
%
%Función Potencia para el contraste unilateral
vmu1=mu:0.01:8; vmu1=vmu1';
vpotencia=zeros(length(vmu1),1);
for i=1:length(vmu1)
    vpotencia(i)=1-normcdf((RC_uni-vmu1(i))/(sigma/(T(1)^0.5)),0,1);
end
figure;
plot(vmu1,vpotencia);
title('Curva de potencia; contraste unilateral');
