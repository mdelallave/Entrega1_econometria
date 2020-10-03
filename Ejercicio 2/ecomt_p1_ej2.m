%% Limpiamos
close('all')
clearvars
clc

%% Import Data
data = xlsread('datosescogidos.xlsx','B4:B6713');
full_eur_usd = data(:,1); % Tipo de cambio euro-us
eur_usd = full_eur_usd(1:end-100);
%% Plots on EUR-USD
figure(1) % Data plot
plot(eur_usd)
title('EUR-USD')

% Mean-std plot
part = 100;
[med,dtip] = graf_m_std(eur_usd,part);

figure(3)
subplot(2,1,1) % ACF
autocorr(eur_usd,20);
title('EUR\_USD ACF')

subplot(2,1,2) % PACF
parcorr(eur_usd,20);
title('EUR\_USD PACF')
pause

%% Taking logarithms
log_eur_usd = 100*log(eur_usd);

figure(4) % Data plot
plot(log_eur_usd)
title('log(EUR-USD)')

figure(5)
subplot(2,1,1) % ACF
autocorr(log_eur_usd,20);
title('EUR\_USD ACF')

subplot(2,1,2) % PACF
parcorr(log_eur_usd,20);
title('EUR\_USD PACF')

pause

%% First logarithmic differences
dlog_eur_usd = log_eur_usd(2:end) - log_eur_usd(1:end-1);

figure(6) % Data plot
plot(dlog_eur_usd)
title('EUR-USD Yield')

figure(7)
subplot(2,1,1) % ACF
autocorr(dlog_eur_usd,20);
title('EUR\_USD Yield ACF')

subplot(2,1,2) % PACF
parcorr(dlog_eur_usd,20);
title('EUR\_USD Yield PACF')

pause

%% Modelization: ARIMA(1,1,1)
model_1 = arima(1,1,1);
estim_model_1 = estimate(model_1,log_eur_usd);
%                   Value      StandardError    TStatistic    PValue 
%                 _________    _____________    __________    _______
% 
%     Constant    0.0060161      0.0065739        0.91515     0.36011
%     AR{1}         0.16646        0.26921        0.61834     0.53635
%     MA{1}        -0.13135        0.27067       -0.48529     0.62747
%     Variance      0.34172      0.0034675         98.549           0

[resid,var_resid,logL] = infer(estim_model_1,log_eur_usd);
std_resid = resid/sqrt(estim_model_1.Variance);
[h,p,jbstat,jbcritval] = jbtest(std_resid,0.01)
figure(8)
subplot(2,2,1);
plot(std_resid);
title('Standardized Residuals')
subplot(2,2,2); % Histogram
histogram(std_resid,10);
title('Standardized Residuals')
subplot(2,2,3); % ACF
autocorr(std_resid,20);
subplot(2,2,4); % PACF
parcorr(std_resid,20);
pause

%% Modelization: ARIMA(1,1,1) x ARIMA(15,0,15)
model_2 = arima('ARlags',[1,15,20],'MAlags',[1,15,20],'D',1);
estim_model_2 = estimate(model_2,log_eur_usd);
%                   Value      StandardError    TStatistic      PValue  
%                 _________    _____________    __________    __________
% 
%     Constant    0.0052939      0.0058095       0.91124         0.36217
%     AR{1}        -0.18481        0.07537       -2.4521        0.014204
%     AR{15}        -0.1403       0.063665       -2.2038        0.027538
%     AR{20}        0.61559       0.069504        8.8569      8.2261e-19
%     MA{1}         0.21023       0.075098        2.7995       0.0051189
%     MA{15}        0.17778       0.064687        2.7483       0.0059899
%     MA{20}        -0.5931        0.06925       -8.5646      1.0844e-17
%     Variance      0.34046      0.0034486        98.725               0

[resid_2,var_resid_2,logL_2] = infer(estim_model_2,log_eur_usd);
std_resid_2 = resid_2/sqrt(estim_model_2.Variance);
[h_2,p_2,jbstat_2,jbcritval_2] = jbtest(std_resid_2,0.01)

figure(9)
subplot(2,2,1);
plot(std_resid_2);
title('Standardized Residuals')
subplot(2,2,2); % Histogram
histogram(std_resid_2,10);
title('Standardized Residuals')
subplot(2,2,3); % ACF
autocorr(std_resid_2,20);
subplot(2,2,4); % PACF
parcorr(std_resid_2,20);
pause

%% Forecast
N = 100;
[log_eur_usd_f,var_err_f] = forecast(estim_model_2,N,'Y0',log_eur_usd);

lower_band = log_eur_usd_f - 1.96*sqrt(var_err_f);
upper_band = log_eur_usd_f + 1.96*sqrt(var_err_f);

figure(10)
plot(6000:6610,log_eur_usd(6000:6610),'Color',[.7,.7,.7],'LineWidth',1.1);
hold on
h1 = plot(6611:6710,lower_band,'r:','LineWidth',2);
plot(6611:6710,upper_band,'r:','LineWidth',2)
h2 = plot(6611:6710,log_eur_usd_f,'k','LineWidth',2);
legend([h1 h2],'Intervalo al 95% de confianza','Predicción puntual',...
	     'Location','NorthWest')
title('log(EUR-USD) Forecast')
hold off
pause
%
psi=impulse(estim_model_2,20);
figure;
subplot(2,1,1)
plot(psi)
title('Impulse response function for log(EUR-USD)')
subplot(2,1,2)
plot(psi(2:end)-psi(1:end-1))
title('Impulse response function for dlog(EUR-USD)')

%% Forecast evaluation
log_full_eur_usd = log(full_eur_usd);
MSE = sum((log_eur_usd_f - log_full_eur_usd(end-99:end)).^2)/N;
RMSE = sqrt(MSE);
MAE = sum(abs(log_eur_usd_f - log_full_eur_usd(end-99:end)));
pMSE = sum(((log_eur_usd_f - log_full_eur_usd(end-99:end))./...
                                   log_full_eur_usd(end-99:end)).^2)/N;
pRMSE = sqrt(pMSE);
pMAE = sum(abs((log_eur_usd_f - log_full_eur_usd(end-99:end))./...
                                log_full_eur_usd(end-99:end)));









