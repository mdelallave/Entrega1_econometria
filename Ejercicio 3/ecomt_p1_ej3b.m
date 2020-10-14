
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
xlabel('Períodos')
ylabel('Cotización')
title('IBEX35');
figure(2);
subplot(2,1,1);
autocorr(ibex,20);
title('ACF del IBEX');
subplot(2,1,2);
parcorr(ibex,20); 
title('PACF del  IBEX');

%% Gráfico media-varianza
part=40;
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

%% Modelo Definitivo
modelo4 = arima('Constant',0,'ARLags',[2,4,5,10],'MALags',[5,10],'D',1);
estmodelo4=estimate(modelo4,libex); 
% Todos los parámetros ahora vuelves a ser significativos al 5%.

[res4,varres4,logL4] = infer(estmodelo4,libex); %inferencia sobre el modelo.
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
% que es significativa

%% Variable impulso
dummyImp = zeros(1569,1);
dummyImp(958) = 1;

%% Modelo Definitivo con impulso
estmodelo_imp = estimate(modelo4,libex,'X', dummyImp); 

[res_imp,varres_imp,logL_imp] = infer(estmodelo_imp,libex); %inferencia sobre el modelo.
resstd_imp = res_imp/sqrt(estmodelo_imp.Variance);
figure;
subplot(1,2,1);
autocorr(resstd_imp,20);
subplot(1,2,2);
parcorr(resstd_imp,20);

% El coeficiente asociado a la variable instrumental NO es significativo, con
% lo que nos quedamos con eL modelo anterior (igual pero excluyendo la
% dummy)

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

plot(1:56+10,f_r(1:50,:))
title('Rendimientos del IBEX simulados para los últimos 56 días del año')
xlabel('t')
ylabel('Rendimiento estimado')

r_acum = sum(f_r,2);

r_acum_pos = r_acum > 0;

n_positivos = sum(r_acum_pos);

ptje_pos = n_positivos/M;

