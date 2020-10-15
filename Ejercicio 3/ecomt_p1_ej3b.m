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
% Nos encontramos que el lag 2 y 10 muestran estructura, ya que se encuentra
% fuera de las bandas.
%% Modelo 3
modelo3=arima('Constant',0,'ARLags',[2,4,5,10],'MALags',[2,4,5,10,12],'D',1);
estmodelo3=estimate(modelo3,libex); 

%     ARIMA(10,1,12) Model (Gaussian Distribution):
%  
%                   Value      StandardError    TStatistic      PValue  
%                 _________    _____________    __________    __________
% 
%     Constant            0             0            NaN             NaN
%     AR{2}        -0.16537      0.037335        -4.4292      9.4584e-06
%     AR{4}       -0.083859      0.031568        -2.6565       0.0078963
%     AR{5}         -0.3431       0.12648        -2.7127       0.0066745
%     AR{10}        0.50662       0.11874         4.2666      1.9845e-05
%     MA{2}         0.12919      0.043082         2.9987       0.0027114
%     MA{4}        0.030318      0.034955        0.86735         0.38575
%     MA{5}         0.28807       0.12177         2.3657        0.017998
%     MA{10}       -0.56656       0.11272        -5.0261      5.0064e-07
%     MA{12}      -0.025487      0.022397         -1.138         0.25514
%     Variance       1.4089      0.026061         54.064               0
%
% MA(4) y el MA(12) no son significativos al 5%.
%% Modelo Definitivo
modelo4 = arima('Constant',0,'ARLags',[2,4,5,10],'MALags',[5,10],'D',1);
estmodelo4=estimate(modelo4,libex); 

% ARIMA(10,1,10) Model (Gaussian Distribution):
%  
%                   Value      StandardError    TStatistic      PValue  
%                 _________    _____________    __________    __________
% 
%     Constant            0             0            NaN             NaN
%     AR{2}        -0.06927      0.015268        -4.5369      5.7087e-06
%     AR{4}       -0.047275      0.015109         -3.129       0.0017539
%     AR{5}        -0.30971       0.12732        -2.4325        0.014996
%     AR{10}        0.50462       0.12175         4.1447       3.402e-05
%     MA{5}         0.24658       0.12197         2.0217        0.043208
%     MA{10}       -0.55782       0.11364        -4.9086      9.1715e-07
%     Variance       1.4166      0.025066         56.516               0
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
dummyImp = zeros(1569,1);
dummyImp(958) = 1;

%% Modelo Definitivo con impulso
estmodelo_imp = estimate(modelo4,libex,'X', dummyImp); 

%   ARIMAX(10,1,10) Model (Gaussian Distribution):
%  
%                   Value      StandardError    TStatistic      PValue  
%                 _________    _____________    __________    __________
% 
%     Constant            0             0            NaN             NaN
%     AR{2}       -0.068329      0.015445        -4.4239      9.6932e-06
%     AR{4}       -0.046697      0.015511        -3.0105       0.0026079
%     AR{5}        -0.33073        0.1315         -2.515        0.011902
%     AR{10}        0.47496       0.12715         3.7354      0.00018744
%     MA{5}          0.2651       0.12624            2.1        0.035729
%     MA{10}       -0.52913       0.11891        -4.4496      8.6011e-06
%     Beta(1)       0.94662       0.76366         1.2396         0.21513
%     Variance       1.4155      0.025156         56.269               0
%
[res_imp,varres_imp,logL_imp] = infer(estmodelo_imp,libex); %inferencia sobre el modelo.
resstd_imp = res_imp/sqrt(estmodelo_imp.Variance);
figure;
subplot(1,2,1);
autocorr(resstd_imp,20);
subplot(1,2,2);
parcorr(resstd_imp,20);

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

