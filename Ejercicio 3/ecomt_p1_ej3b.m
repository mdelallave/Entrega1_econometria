
%% Limpiamos
close('all')
clear all
clc

%% Importar datos
data = xlsread('ibex.xlsx','B3:C1558');
ibex = data; 

%% Hacemos un grafico con Ibex
figure(1);
plot(ibex);
title('Ibex');
figure(2);
subplot(2,1,1);
autocorr(ibex,20);
title('ACF del IBEX');
subplot(2,1,2);
parcorr(ibex,20); 
title('PACF del  IBEX');

%% Gráfico media-varianza
part=20;
[med,dtip] = graf_m_std(ibex,part);

%% Logaritmos
% Transformación logarítmica (rendimientos)
libex = log(ibex);

% Estudiamos su estacionariedad
figure(10);
subplot(2,1,1);
autocorr(libex,20);
title('ACF del log del IBEX');
subplot(2,1,2);
parcorr(libex,20); 
title('PACF del log del IBEX');

% Vemos como la ACF toma valores muy altos que no tienden a desaparecer.
% Indica ausencia de estacionariedad.

%% Tomamos primeras diferencias log y analizamos estacionariedad.
dlibex = libex(2:end) - libex(1:end-1);
figure(11);
plot(dlibex);
title('Rendimiento del Ibex');
figure(4);
subplot(2,1,1);
autocorr(dlibex,20); 
title('ACF del Rendimiento del Ibex');
subplot(2,1,2);
parcorr(dlibex,20);
title('PACF del Rendimiento del Ibex');

% LA ACF ya no presenta infinitos valores significativos. La serie presenta
% estacionariedad. Sin embargo, dlibex presenta un valor muy grande en la
% observación 958 que mas tarde trataremos de controlar mediante una variable dummy.

% %% Modelo 1
% l_model1=arima('ARLags',[4 5 12],'MALags',[4 5 12],'D',1);
% l_estmodelo1=estimate(l_model1,libex); 
% %     ARIMA(12,1,12) Model (Gaussian Distribution):
% %  
% %                   Value      StandardError    TStatistic      PValue  
% %                 _________    _____________    __________    __________
% % 
% %     Constant    0.0086269      0.054755        0.15755         0.87481
% %     AR{4}        -0.32981       0.12231        -2.6965       0.0070075
% %     AR{5}        -0.59261       0.10378          -5.71      1.1295e-08
% %     AR{12}       0.074844      0.071569         1.0458         0.29567
% %     MA{4}         0.30966       0.12364         2.5047        0.012257
% %     MA{5}         0.56652       0.10675         5.3072      1.1133e-07
% %     MA{12}       -0.14229      0.072973        -1.9499        0.051194
% %     Variance       1.4218      0.028957           49.1               0
% 
% %% Hacemos inferencia
% [res,varres,logL] = infer(l_estmodelo1,libex); %inferencia sobre el modelo.
% l_resstd1 = res/sqrt(l_estmodelo1.Variance);
% %test jarque bera. h=1 rechazo normalidad.
% figure;
% subplot(2,2,1);
% plot(l_resstd1); 
% title('Residuos estandarizados modelo 1');
% subplot(2,2,2);
% histogram(l_resstd1,10);
% title('Residuos Estandarizados modelo 1')
% subplot(2,2,3);
% autocorr(l_resstd1,20);
% subplot(2,2,4);
% parcorr(l_resstd1,20);
% 
% %% Eliminamos la constante y los retardos no significativos
% l_model2=arima('Constant',0,'ARLags',[4 5],'MALags',[4 5],'D',1);
% l_estmodelo2=estimate(l_model2,libex);
% %     ARIMA(5,1,5) Model (Gaussian Distribution):
% %  
% %                   Value       StandardError    TStatistic      PValue  
% %                 __________    _____________    __________    __________
% % 
% %     Constant             0             0             NaN            NaN
% %     AR{4}          -0.7013      0.075381         -9.3034     1.3606e-20
% %     AR{5}        -0.041545      0.087844        -0.47294        0.63626
% %     MA{4}          0.64626      0.081449          7.9346     2.1121e-15
% %     MA{5}       -0.0099157      0.095269        -0.10408        0.91711
% %     Variance        1.4306      0.027522          51.979              0
% 
% %% Hacemos inferencia
% [res,varres,logL] = infer(l_estmodelo2,libex); %inferencia sobre el modelo.
% l_resstd2 = res/sqrt(l_estmodelo2.Variance);
% [h,p,jbstat,critval] = jbtest(l_resstd2, 0.01);
% %test jarque bera. h=1 rechazo normalidad.
% figure;
% subplot(2,2,1);
% plot(l_resstd2); 
% title('Residuos estandarizados modelo 1');
% subplot(2,2,2);
% histogram(l_resstd2,10);
% title('Residuos Estandarizados modelo 1')
% subplot(2,2,3);
% autocorr(l_resstd2,20);
% subplot(2,2,4);
% parcorr(l_resstd2,20);
% 
% %% Reajustamos de nuevo el modelo
% l_model3=arima('Constant',0,'ARLags',[4 5 10 11 12],'MALags',[5 11 10],'D',1);
% l_estmodelo3=estimate(l_model3,libex);
% %     ARIMA(12,1,11) Model (Gaussian Distribution):
% %  
% %                   Value      StandardError    TStatistic      PValue  
% %                 _________    _____________    __________    __________
% % 
% %     Constant            0             0            NaN             NaN
% %     AR{4}       -0.034189      0.014759        -2.3164        0.020534
% %     AR{5}        -0.50935       0.12783        -3.9845      6.7612e-05
% %     AR{10}        0.29324       0.13601          2.156        0.031084
% %     AR{11}       0.019387      0.048533        0.39947         0.68955
% %     AR{12}      -0.062996      0.015584        -4.0423      5.2927e-05
% %     MA{5}         0.43789        0.1251         3.5004       0.0004646
% %     MA{10}       -0.36482        0.1305        -2.7955       0.0051823
% %     MA{11}      -0.061438      0.053563         -1.147         0.25137
% %     Variance       1.4122      0.028051         50.342               0
% %% Hacemos inferencia
% [res,varres,logL] = infer(l_estmodelo3,libex); %inferencia sobre el modelo.
% l_resstd3 = res/sqrt(l_estmodelo3.Variance);
% [h,p,jbstat,critval] = jbtest(l_resstd3, 0.01);
% %test jarque bera. h=1 rechazo normalidad.
% figure;
% subplot(2,2,1);
% plot(l_resstd3); 
% title('Residuos estandarizados modelo 1');
% subplot(2,2,2);
% histogram(l_resstd3,10);
% title('Residuos Estandarizados modelo 1')
% subplot(2,2,3);
% autocorr(l_resstd3,20);
% subplot(2,2,4);
% parcorr(l_resstd3,20);

%% Último modelo de Diego
modelo_Diego = arima('Constant',0,'ARLags',[4 5 10 12],'MALags',[5 10],'D',1);
estmodelo_Diego = estimate(modelo_Diego,libex);
%     ARIMA(12,1,10) Model (Gaussian Distribution):
%  
%                   Value      StandardError    TStatistic      PValue  
%                 _________    _____________    __________    __________
% 
%     Constant            0             0            NaN             NaN
%     AR{4}       -0.043813      0.015318        -2.8602       0.0042334
%     AR{5}        -0.44163        0.1362        -3.2425       0.0011848
%     AR{10}        0.35273        0.1389         2.5395        0.011102
%     AR{12}      -0.064831      0.016202        -4.0015      6.2929e-05
%     MA{5}         0.36578       0.13268         2.7569        0.005836
%     MA{10}       -0.41847       0.13239        -3.1609       0.0015728
%     Variance       1.4175      0.026045         54.427               0
%% Hacemos inferencia sobre el modelo de Diego
[res,varres,logL_imp] = infer(estmodelo_Diego,libex); %inferencia sobre el modelo.
l_resstd4 = res/sqrt(estmodelo_Diego.Variance);
[AIC_Diego,BIC_Diego] = aicbic(logL_imp,6,1556)
figure(2);
subplot(2,2,1);
plot(l_resstd4); 
title('Residuos estandarizados modelo 1');
subplot(2,2,2);
histogram(l_resstd4,10);
title('Residuos Estandarizados modelo 1')
subplot(2,2,3);
autocorr(l_resstd4,20);
subplot(2,2,4);
parcorr(l_resstd4,20);

%% Modelo 4 Carlos (DEFINITIVO)
modelo4 = arima('Constant',0,'ARLags',[2,4,5,10],'MALags',[2,5,10],'D',1);
estmodelo4=estimate(modelo4,libex); 
% Todos los parámetros ahora vuelves a ser significativos al 5%.

[res4,varres4,logL4] = infer(estmodelo4,libex); %inferencia sobre el modelo.
resstd4 = res4/sqrt(estmodelo4.Variance);
[AIC_Carlos,BIC_Carlos] = aicbic(logL4,9,1556);
figure;
subplot(2,2,1);
plot(resstd4); 
title('Residuos estandarizados modelo 4');
subplot(2,2,2);
histogram(resstd4,10);
title('Residuos Estandarizados modelo 4')
subplot(2,2,3);
autocorr(resstd4,20);
subplot(2,2,4);
parcorr(resstd4,20);

% Este modelo es el mejor y es el que usamos para la predicción mediante
% simulación. No obstante, le añadimos la variable dummy a continuación,
% que es significativa

%% Función impulso
dummyImp = zeros(1569,1);
dummyImp(958) = 1;

%% Modelo Definitivo con impulso
estmodelo_imp = estimate(modelo4,libex,'X', dummyImp); 

% El coeficiente asociado a la variable instrumental es significativo, con
% lo que nos quedamos con este modelo

%% Predicciones (simulacion alternativa, no funciona)
var_resid = estmodelo4.Variance;

M = 1.0e4; % número de trayectorias simuladas
n = 56; % número de observaciones por trayectoria
epsilon = normrnd(0,sqrt(var_resid),M,n + 12); 

% Parameters
phi_2 = estmodelo4.AR{2};
phi_4 = estmodelo4.AR{4};
phi_5 = estmodelo4.AR{5};
phi_10= estmodelo4.AR{10};
theta_2 = estmodelo4.MA{2};
theta_5 = estmodelo4.MA{5};
theta_10 = estmodelo4.MA{10};

f_r = zeros(M,n + 12);

for i = 1:M
    epsilon(i,1:12) = res4(end-11:end);
    f_r(i,1:12) = dlibex(end-11:end); % Primeras diferencias logarítmicas
    for j = 13:n + 12
    f_r(i,j) = phi_2*f_r(i,j - 2) + phi_4*f_r(i,j - 4) + ...
        phi_5*f_r(i,j - 5) + phi_10*f_r(i,j - 10) - ... 
        theta_2*epsilon(i,j - 2) - theta_5*epsilon(i,j - 5) - ...
        theta_10*epsilon(i,j - 10) + epsilon(i,j); 
    end
end

media_r = mean(f_r);
var_r = var(f_r);

