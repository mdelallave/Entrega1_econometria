%% Limpiamos
close('all')
clear all
clc

%% Importar datos
data = xlsread('ibex.xlsx','B3:C1558');
ibex = data; 

%% Hacemos un gráfico con Ibex
figure(1);
plot(ibex);
title('Ibex');
figure(2)
subplot(2,1,1);
autocorr(ibex,20);
title('ACF del IBEX');
subplot(2,1,2);
parcorr(ibex,20); 
title('PACF del  IBEX');

%% Analisis de estacionariedad.
% Aplicamos un gráfico media-varianza
part=100;
[med,dtip] = graf_m_std(ibex,part);
% No queda claro si los datos en nivel presentan estacionariedad usando
% el gráfico media varianza. Sin embargo, apicamos logaritmos por sus
% buenas propiedades estadísticas.

%% Manipulación de los datos
% Transformación logarítmica (rendimientos)
libex = 100*log(ibex);
figure(3)
plot(libex)

% Estudiamos su estacionariedad
figure(4);
subplot(2,1,1);
autocorr(libex,20);
title('ACF del log del IBEX');
subplot(2,1,2);
parcorr(libex,20); 
title('PACF del log del IBEX');

% Vemos como la ACF toma valores muy altos que no tienden a desaparecer.
% Indica ausencia de estacionariedad.

%% Tomamos primeras diferencias y analizamos estacionariedad.
dlibex=libex(2:end)-libex(1:end-1);
figure(5);
plot(dlibex);
title('Rendimiento logarítmico del Ibex');
figure(6);
subplot(2,1,1);
autocorr(dlibex,20); 
title('ACF del Rendimiento logarítmico del Ibex');
subplot(2,1,2);
parcorr(dlibex,20);
title('PACF del Rendimiento logarítmico del Ibex');

% LA ACF ya no presenta infinitos valores significativos. La serie presenta
% estacionariedad. Sin embargo, dlibex presenta un valor muy grande en la
% observación 958 que mas tarde trataremos de controlar mediante una variable dummy.

%% Modelo 1
modelo1=arima('ARlags',[4,5,12], 'MAlags', [4,5,12], 'D', 1);
estmodelo1=estimate(modelo1,libex); 

%     ARIMA(12,1,12) Model (Gaussian Distribution):
%  
%                   Value      StandardError    TStatistic      PValue  
%                 _________    _____________    __________    __________
% 
%     Constant    0.0086261      0.054755        0.15754         0.87482
%     AR{4}        -0.32981       0.12231        -2.6965       0.0070076
%     AR{5}        -0.59261       0.10378        -5.7101      1.1291e-08
%     AR{12}       0.074843      0.071569         1.0457         0.29568
%     MA{4}         0.30966       0.12363         2.5046        0.012257
%     MA{5}         0.56652       0.10674         5.3073      1.1129e-07
%     MA{12}       -0.14228      0.072972        -1.9498        0.051195
%     Variance       1.4218      0.028957           49.1               0
%
% Nos encontramos que la constante, el AR(12) no son significativas al 5%.
% Sin embargo, el MA(12) está en los límites de significtividad. Lo
% dejamos.
%% Modelo 2
modelo2=arima('Constant',0,'ARLags',[4,5],'MALags',[4,5,12],'D',1);
estmodelo2=estimate(modelo2,libex); 

%   ARIMA(5,1,12) Model (Gaussian Distribution):
%  
%                   Value      StandardError    TStatistic      PValue  
%                 _________    _____________    __________    __________
% 
%     Constant            0             0            NaN             NaN
%     AR{4}        -0.22427      0.070015        -3.2031       0.0013596
%     AR{5}        -0.66355      0.064308        -10.318      5.8253e-25
%     MA{4}         0.19703      0.072331          2.724       0.0064488
%     MA{5}         0.63405      0.070089         9.0463      1.4783e-19
%     MA{12}      -0.062362       0.01805        -3.4549       0.0005504
%     Variance       1.4227       0.02798         50.845               0
%
% Todo es significativo al 5%.

%% Hacemos inferencia
[res2,varres2,logL2] = infer(estmodelo2,libex); %inferencia sobre el modelo.
resstd2 = res2/sqrt(estmodelo2.Variance);
[h,p,jbstat,critval] = jbtest(resstd2, 0.01) 
%test jarque bera. h=1 rechazo normalidad.
figure;
subplot(2,2,1);
plot(resstd2); 
title('Residuos estandarizados modelo 1');
subplot(2,2,2);
histogram(resstd2,10);
title('Residuos Estandarizados modelo 1')
subplot(2,2,3);
autocorr(resstd2,16);
subplot(2,2,4);
parcorr(resstd2,16);
% Nos encontramos que el lag 10 muestran estructura, ya que se encuentra
% fuera de las bandas. 
%% Modelo 3
modelo3=arima('Constant',0,'ARLags',[4,5,10],'MALags',[4,5,10,12],'D',1);
estmodelo3=estimate(modelo3,libex); 

%  ARIMA(10,1,12) Model (Gaussian Distribution):
%  
%                   Value      StandardError    TStatistic      PValue  
%                 _________    _____________    __________    __________
% 
%     Constant            0             0            NaN             NaN
%     AR{4}        -0.20214      0.066731        -3.0292       0.0024517
%     AR{5}        -0.49088       0.14727        -3.3333      0.00085827
%     AR{10}        0.21739       0.13828         1.5722         0.11591
%     MA{4}         0.17407      0.068263           2.55        0.010773
%     MA{5}         0.41753       0.14688         2.8426       0.0044742
%     MA{10}       -0.28372       0.13335        -2.1277        0.033364
%     MA{12}      -0.068092      0.018238        -3.7335      0.00018886
%     Variance       1.4159      0.027992         50.582               0
%
% AR(10) no son significativos al 5%.
%% Modelo Definitivo
modelo4 = arima('Constant',0,'ARLags',[4,5],'MALags',[4,5,10,12],'D',1);
estmodelo4=estimate(modelo4,libex); 

%  ARIMA(5,1,12) Model (Gaussian Distribution):
%  
%                  Value      StandardError    TStatistic      PValue  
%                 ________    _____________    __________    __________
% 
%     Constant           0             0            NaN             NaN
%     AR{4}       -0.18211      0.058309        -3.1232       0.0017893
%     AR{5}       -0.71431      0.053769        -13.285      2.8405e-40
%     MA{4}        0.15876      0.060144         2.6397       0.0082989
%     MA{5}         0.6403       0.05951          10.76      5.3401e-27
%     MA{10}      -0.06848      0.023412         -2.925       0.0034449
%     MA{12}      -0.06214      0.017157        -3.6219      0.00029242
%     Variance      1.4167      0.027528         51.466               0
%
% Todos los parámetros ahora vuelves a ser significativos al 5%.
%% Hacemos inferencia con el modelo 4
[res4,varres4,logL4] = infer(estmodelo4,libex);%inferencia sobre el modelo.
resstd4 = res4/sqrt(estmodelo4.Variance);
figure;
subplot(2,2,1);
plot(resstd4); 
title('Residuos estandarizados modelo 4');
subplot(1,2,1);
autocorr(resstd4,20);
subplot(1,2,2);
parcorr(resstd4,20);

% Este modelo es el mejor y es el que usamos para la predicción mediante
% simulación. No obstante, le añadimos la variable dummy a continuación,
% que no es significativa

%% Variable impulso
dummyImp = zeros(1574,1);
dummyImp(958) = 1;

%% Modelo Definitivo con impulso
estmodelo_imp = estimate(modelo4,libex,'X', dummyImp); 

% ARIMAX(5,1,12) Model (Gaussian Distribution):
%  
%                   Value      StandardError    TStatistic      PValue  
%                 _________    _____________    __________    __________
% 
%     Constant            0             0            NaN             NaN
%     AR{4}        -0.39601      0.090914        -4.3559      1.3255e-05
%     AR{5}         0.36551      0.089973         4.0624      4.8572e-05
%     MA{4}         0.35387       0.08632         4.0996      4.1391e-05
%     MA{5}        -0.45351      0.093632        -4.8435      1.2755e-06
%     MA{10}       0.010118      0.021375        0.47334         0.63597
%     MA{12}      -0.043535      0.019912        -2.1864        0.028788
%     Beta(1)        1.8546       0.52268         3.5482       0.0003879
% El MA(10) no es significativo al 5%. (Al hacer inferencia con el modelo
% resultante, el lag 17 muestra estructura)
%% Modelo 
modelo5 = arima('Constant',0,'ARLags',[4,5,17],'MALags',[4,5,12,17],'D',1);
estmodelo5=estimate(modelo5,libex); 
% Modelo optimo
estmodelo_imp2 = estimate(modelo5,libex,'X', dummyImp); 

% Inferencia con el modelo optimo
[res_imp2,varres_imp2,logL_imp2] = infer(estmodelo_imp2,libex); %inferencia sobre el modelo.
resstd_imp2 = res_imp2/sqrt(estmodelo_imp2.Variance);
figure;
subplot(1,2,1);
autocorr(resstd_imp2,20);
subplot(1,2,2);
parcorr(resstd_imp2,20);

% Podemos comprobar que el  coeficiente asociado a la variable instrumental
% NO es significativo, con lo que nos quedamos con eL modelo anterior
%(igual pero excluyendo la dummy)

%% Predicciones 
var_resid = estmodelo4.Variance;

M = 1.0e4; % número de trayectorias simuladas
n = 56; % número de observaciones por trayectoria
epsilon = normrnd(0,sqrt(var_resid),M,n + 10); 

% Parameters
phi_2 = estmodelo4.AR{2};
phi_4 = estmodelo4.AR{4};
phi_5 = estmodelo4.AR{5};
phi_10= estmodelo4.AR{10};
theta_5 = estmodelo4.MA{5};
theta_10 = estmodelo4.MA{10};

f_r = zeros(M,n + 10);

for i = 1:M
    epsilon(i,1:10) = res4(end-9:end);
    f_r(i,1:10) = dlibex(end-9:end); % Primeras diferencias logarítmicas
    for j = 11:n + 10
    f_r(i,j) = phi_2*f_r(i,j - 2) + phi_4*f_r(i,j - 4) + ...
        phi_5*f_r(i,j - 5) + phi_10*f_r(i,j - 10) - ... 
        theta_5*epsilon(i,j - 5) - theta_10*epsilon(i,j - 10) + epsilon(i,j); 
    end
end
% Graficamos las simulaciones y su media.
figure;
plot(1:n+10,f_r(1:50,:))
hold on
plot(1:n+10,mean(f_r),'k','LineWidth',1.2)
legend('50 trayectorias simuladas','media (línea negra)')
xlabel('t')
ylabel('Rendimientos del IBEX35')
title('Rendimientos del IBEX simulados para los últimos 56 días de negociación')
hold off

r_acum = sum(f_r,2);

r_acum_pos = r_acum > 0;

n_positivos = sum(r_acum_pos);

ptje_pos = n_positivos/M;

