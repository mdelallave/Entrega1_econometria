%% Limpiamos
close('all')
clearvars
clc

%% Importar datos
% data = csvread('GDPC1.csv');
% GPD = data; 
% En el mac no me va la función de importar datos con xlsread. Y con csvread no sé hacerlo. Probad si os
% funciona.
%% Hacemos un gráfico con GPD
figure(1);
plot(gdpc);
title('gdpc');
%% Particionamos la muestra en 100 submuestras, 
% y hacemos el grafico media desviacion tipica
part=100;
[med,dtip] = graf_m_std(gdpc,part);
% Los datos en nivel no presentan estacionariedad

%% Manipulación de los datos
% Transformación logarítmica (rendimientos)
lgdpc = 100*log(gdpc);

% Estudiamos su estacionariedad
figure(2);
subplot(2,1,1);
autocorr(lgdpc,20);
title('ACF del log del GDPC');
subplot(2,1,2);
parcorr(lgdpc,20); 
title('PACF del log del GDPC');

% Vemos como la ACF toma valores muy altos que no tienden a desaparecer.
% Indica ausencia de estacionariedad.

%% Tomamos primeras diferencias
dlgdpc=lgdpc(2:end)-lgdpc(1:end-1);
figure(3);
plot(dlgdpc);
title('Tasa de variación de GPDC');
figure(4);
subplot(2,1,1);
autocorr(dlgdpc,20); 
title('ACF de la tasa de variación de GPDC');
subplot(2,1,2);
parcorr(dlgdpc,20);
title('PACF de la tasa de variación de GPDC');
% LA ACF ya no presenta valores positivos. La serie ya presenta
% estacionariedad.

%% Modelo
% Empezamos con un modelo ARIMA(1,1,1)
modelo=arima(1,1,1);
estmodelo=estimate(modelo,dlgdpc); 
% estimamos el modelo.La constante no es significativa.

%% Hacemos inferencia
[res,varres,logL] = infer(estmodelo,dlgdpc); %inferencia sobre el modelo.
resstd = res/sqrt(estmodelo.Variance);
[h,p,jbstat,critval]=jbtest(resstd,0.01) 
%test jarque bera. h =0 -> no puedo rechazar h0 de normalidad. h=1 rechazo normalidad.
% rechazo normalidad.
figure;
subplot(2,2,1);
plot(resstd); 
title('Residuos estandarizados');
subplot(2,2,2);
histogram(resstd,10);
title('Residuos Estandarizados')
subplot(2,2,3);
autocorr(resstd,16);
subplot(2,2,4);
parcorr(resstd,16);
%% Modelo 2
modelo2=arima( 'ARLags',[1,2,9,10,11,15],'MALags',[1,13],'D',1);
estmodelo2=estimate(modelo2,dlgdpc); 
% estimamos el modelo.La constante no es significativa.

%% Hacemos inferencia
[res2,varres2,logL2] = infer(estmodelo2,dlgdpc); %inferencia sobre el modelo.
resstd2 = res2/sqrt(estmodelo2.Variance);
[h2,p2,jbstat2,critval2]=jbtest(resstd2,0.01) 
%test jarque bera. h =0 -> no puedo rechazar h0 de normalidad. h=1 rechazo normalidad.
% rechazo normalidad.
figure;
subplot(2,2,1);
plot(resstd2); 
title('Residuos estandarizados modelo 2');
subplot(2,2,2);
histogram(resstd2,10);
title('Residuos Estandarizados modelo 2')
subplot(2,2,3);
autocorr(resstd2,16);
subplot(2,2,4);
parcorr(resstd2,16);

% ¿Está el modelo bien especificado? No encuentro de reducir la
% autocorrelación en los lags cercanos a la barrera.
%% Hacemos una predicción
