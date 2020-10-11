%% Limpiamos
close('all')
clearvars
clc

%% Importar datos
data = xlsread('data.csv','B4:D6613');
eurous = data(:,1); % Tipo de cambio euro-us
repsol = data(:,2); % Repsol
crude = data(:,3);  % Futuros del crudo

%% Parámetros iniciales
N = 100; % Cantidad de datos omitidos.
alpha = 0.05; % Nivel de significación.
Z = norminv(alpha/2,0,1);
a = 100-alpha*100; a = string(a); % Alpha como caracter.
T = length(data(1:end-N,1)); % Cantidad de datos.

%% Manipulación de los datos
% Transformación logarítmica (rendimientos)
leurous = 100*log(eurous(1:end-N)); lrepsol = 100*log(repsol(1:end-N)); 
lcrude = 100*log(crude(1:end-N)); ldata = 100*log(data(1:end-N,:));
% Estudiamos su ACF y PACF
figure(1);
subplot(3,2,1)
autocorr(leurous,20);
title('ACF del log del EURO-US');
subplot(3,2,2)
parcorr(leurous,20);
title('PACF del log del EURO-US');

subplot(3,2,3);
autocorr(lrepsol,20);
title('ACF del log de Repsol');
subplot(3,2,4);
parcorr(lrepsol,20);
title('PACF del log de Repsol');

subplot(3,2,5);
autocorr(lcrude,20);
title('ACF del log del crudo');
subplot(3,2,6);
parcorr(lcrude,20);
title('PACF del log del crudo');
% Necesitamos tomar diferencias: ACF decrece lenta y linealmente

%% Tomamos diferencias
dleurous = leurous(2:end)-leurous(1:end-1);
dlrepsol = lrepsol(2:end)-lrepsol(1:end-1);
dlcrude = lcrude(2:end)-lcrude(1:end-1);
dldata = ldata(2:end,:)-ldata(1:end-1,:);

figure(2);
subplot(3,2,1)
autocorr(dleurous,20);
title('ACF del rendimiento del EURO-US');
subplot(3,2,2)
parcorr(dleurous,20);
title('PACF del rendimiento del EURO-US');

subplot(3,2,3);
autocorr(dlrepsol,20);
title('ACF del rendimiento de Repsol');
subplot(3,2,4);
parcorr(dlrepsol,20);
title('PACF del rendimiento de Repsol');

subplot(3,2,5);
autocorr(dlcrude,20);
title('ACF del rendimiento del crudo');
subplot(3,2,6);
parcorr(dlcrude,20);
title('PACF del rendimiento del crudo');
%% Modelo
% Planteamos inicialmente un modelo ARIMA(1,1,1)

for i = 1:length(ldata(1,:))
    estimacion(i) = estimate(arima(1,1,1),ldata(:,i));
end
% EURO-US nada significativo. Repsol y Crude cte. no significativa.

%% Análisis de los residuos
res = zeros(length(ldata(:,1)),length(ldata(1,:)));
varres = zeros(length(ldata(:,1)),length(ldata(1,:)));
logL = zeros(length(ldata(1,:)),1);
resstd = zeros(length(ldata(:,1)),length(ldata(1,:)));
h = zeros(length(ldata(1,:)),1);
p = zeros(length(ldata(1,:)),1);
jbstat = zeros(length(ldata(1,:)),1);
critval = zeros(length(ldata(1,:)),1);

for i = 1:length(ldata(1,:))
    [res(:,i),varres(:,i),logL(i)] = infer(estimacion(i),ldata(:,i));
    resstd(:,i) = res(:,i)/sqrt(estimacion(i).Variance);
    [h(i),p(i),jbstat(i),critval(i)] = jbtest(resstd(:,i),alpha);
end
h, p, jbstat, critval
% Rechazamos la hipótesis de normalidad para todas las series

%% Gráficos de los residuos
nombre = ["EURO-US","Repsol","crudo"];
for i = 1:length(estimacion)
    figure(2+i);
    subplot(2,2,1);
    plot(resstd(:,i));
    title(['Residuos estandarizados de', nombre(i)]);
    
    subplot(2,2,2);
    histogram(resstd(:,i),20);
    title(['Residuos estandarizados de ', nombre(i)]);
    
    subplot(2,2,3);
    autocorr(resstd(:,i),20);
    
    subplot(2,2,4);
    parcorr(resstd(:,i),20);
end

% Comprobación de los residuos
banda_confianza = 0 - norminv(alpha/2,0,1) * 1/sqrt(T);
significativos_ma = zeros(21,length(estimacion));
lags_significativos_ma = zeros(20,length(estimacion));
significativos_ma = zeros(21,length(estimacion));
lags_significativos_ma = zeros(20,length(estimacion));
IND = 1:20; % Índice

for i = 1:length(estimacion)
    significativos_ma(:,i) = abs(autocorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ma(:,i) = IND'.*significativos_ma(2:end,i);
    significativos_ar(:,i) = abs(parcorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ar(:,i) = IND'.*significativos_ar(2:end,i);
end
lags_significativos_ar
lags_significativos_ma
% EURO-US: AR(6) y MA(6)
% Repsol : AR([2,3,5]) y MA([2,3,5])
% Crude  : AR([1,2,4,8,15]) y MA([1,2,4,8,15])

%% Reajustamos el modelo
% Eliminamos la constante en todos porque sale no significativa.
% EURO-US 
estimacion(1) = estimate(arima('Constant',0,'ARlags',6,'MAlags',6,...
                'D',1),ldata(:,1));
            
% Repsol
estimacion(2) = estimate(arima('Constant',0,'ARlags',[2,3,5],'MAlags',...
                [2,3,5],'D',1),ldata(:,2));

% Crude            
estimacion(3) = estimate(arima('ARlags',[1,2,4,8,15],...
                'MAlags',[1,2,4,8,15],'D',1),ldata(:,3)); 

% EURO-US y Repsol todos parámetros significativos. Crude no
% significativos: AR(15), MA(15) y la cte. pero haciendo pruebas vemos que
% llegamos a un ruido blanco si mantenemos la constante.
%% Volvemos a estudiar los residuos

for i = 1:length(ldata(1,:))
    [res(:,i),varres(:,i),logL(i)] = infer(estimacion(i),ldata(:,i));
    resstd(:,i) = res(:,i)/sqrt(estimacion(i).Variance);
    [h(i),p(i),jbstat(i),critval(i)] = jbtest(resstd(:,i),alpha);
end
h, p, jbstat, critval
% Rechazamos la hipótesis de normalidad para todas las series

%% Gráficas
for i = 1:length(estimacion)
    figure(5+i);
    subplot(2,2,1);
    plot(resstd(:,i));
    title(['Residuos estandarizados de', nombre(i)]);
    
    subplot(2,2,2);
    histogram(resstd(:,i),20);
    title(['Residuos estandarizados de ', nombre(i)]);
    
    subplot(2,2,3);
    autocorr(resstd(:,i),20);
    
    subplot(2,2,4);
    parcorr(resstd(:,i),20);
end

% Comprobación de los residuos
for i = 1:length(estimacion)
    significativos_ma(:,i) = abs(autocorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ma(:,i) = IND'.*significativos_ma(2:end,i);
    significativos_ar(:,i) = abs(parcorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ar(:,i) = IND'.*significativos_ar(2:end,i);
end
lags_significativos_ar
lags_significativos_ma

% EURO-US y Crude ruido blanco. Repsol significativo en AR(1) y MA(1).

%% Reajustamos el modelo

% Repsol
estimacion(2) = estimate(arima('Constant',0,'ARlags',[1,2,3,5],'MAlags',...
                [1,2,3,5],'D',1),ldata(:,2));

% Crude. Eliminamos los lags no significativos            
estimacion(3) = estimate(arima('ARlags',[1,2,4,8],...
                'MAlags',[1,2,4,8],'D',1),ldata(:,3)); 
            
for i = 1:length(estimacion)
    significativos_ma(:,i) = abs(autocorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ma(:,i) = IND'.*significativos_ma(2:end,i);
    significativos_ar(:,i) = abs(parcorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ar(:,i) = IND'.*significativos_ar(2:end,i);
end
lags_significativos_ar
lags_significativos_ma

% Todos los parámetros son significativos excepto Repsol: MA(1) y 
% Crude: AR(2). Respecto a los residuos, Repsol sigue siendo significativo 
% en 1 a pesar de incluir los retardos en el modelo.

%% Reajustamos el modelo

% Crude. Eliminamos los lags no significativos            
estimacion(3) = estimate(arima('ARlags',[1,4,8],...
                'MAlags',[1,2,4,8],'D',1),ldata(:,3)); 

for i = 1:length(estimacion)
    significativos_ma(:,i) = abs(autocorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ma(:,i) = IND'.*significativos_ma(2:end,i);
    significativos_ar(:,i) = abs(parcorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ar(:,i) = IND'.*significativos_ar(2:end,i);
end
lags_significativos_ar
lags_significativos_ma

% Todos los parámetros son significativos excepto la cte. del Crude.
% Haciendo pruebas, vemos que de esta manera los residuos de Crude se
% comportan como ruido blanco, por lo que preferimos dejar la cte.
% Respecto a las correlaciones de los residuos, Repsol sigue siendo 
% significativo en 1.

%% Predicción y bandas de confianza
% Inicialización de variables
ldata_f = zeros(N,length(estimacion)); 
dt_errorf = zeros(N,length(estimacion));
b_inferior = zeros(N,length(estimacion));
b_superior = zeros(N,length(estimacion));


% Predicción estática (a 1 día)
for i = 1:length(estimacion)
    for j = 1:N
        [ldata_f(j,i),dt_errorf(j,i)] = forecast(estimacion(i),1,'Y0',...
            ldata(:,i));
        b_inferior(j,i) = ldata_f(j,i) - Z*sqrt(dt_errorf(j,i));
        b_superior(j,i) = ldata_f(j,i) + Z*sqrt(dt_errorf(j,i));
        ldata = 100*log(data(1:end-N+j,:));
    end
    ldata = 100*log(data(1:end-N,:));
end
% Gráfica
M = round(1.5*N); % Puntos a mostrar
for i = 1:length(estimacion)
    figure(8);
    subplot(3,1,i);
    plot(ldata(end-M:end,i),'Color',[.7,.7,.7]);
    hold on
    h1 = plot(length(ldata(end-M:end,i)):(N-1+length(ldata(end-M:end,i))),...
        b_inferior(:,i),'r:','LineWidth',2);
    plot(length(ldata(end-M:end,i)):(N-1+length(ldata(end-M:end,i))),...
        b_superior(:,i),'r:','LineWidth',2)
    h2 = plot(length(ldata(end-M:end,i)):(N-1+length(ldata(end-M:end,i))),...
        ldata_f(:,i),'k','LineWidth',2);
      legend([h1 h2], a, 'Predicción a 1 día', 'Location','NorthWest')
    title(['Predicción estática del ' nombre(i)])
    hold off
end

%% Predicción dinámica (a k días, con k = 1,2,...,N)
ldata_f2 = zeros(N,length(estimacion)); 
dt_errorf2 = zeros(N,length(estimacion));
for i = 1:length(estimacion)
    [ldata_f2(:,i),dt_errorf2(:,i)] = forecast(estimacion(i),N,'Y0',...
        ldata(:,i));
    b_inferior(:,i) = ldata_f2(:,i) - Z*sqrt(dt_errorf2(:,i));
    b_superior(:,i) = ldata_f2(:,i) + Z*sqrt(dt_errorf2(:,i));
end

for i = 1:length(estimacion)
    figure(9);
    subplot(3,1,i);
    plot(ldata(end-M:end,i),'Color',[.7,.7,.7]);
    hold on
    h1 = plot(length(ldata(end-M:end,i)):(N-1+length(ldata(end-M:end,i))),...
        b_inferior(:,i),'r:','LineWidth',2);
    plot(length(ldata(end-M:end,i)):(N-1+length(ldata(end-M:end,i))),...
        b_superior(:,i),'r:','LineWidth',2)
    h2 = plot(length(ldata(end-M:end,i)):(N-1+length(ldata(end-M:end,i))),...
        ldata_f2(:,i),'k','LineWidth',2);
      legend([h1 h2], a, 'Predicción puntual', 'Location','NorthWest')
    title(['Predicción dinámica del ' nombre(i)])
    hold off
end

%% Comparación valores reales
ldata_all = 100*log(data);
for i = 1:length(estimacion)
    figure(10);
    subplot(3,1,i);
    plot(ldata_all(end-M-N:end,i));
    title(['Valor real del ' nombre(i)])
end

% Error predicción
err_rel = zeros(N,length(estimacion));
for i = 1:length(estimacion)
    err_rel(:,i) = 100*abs(ldata_f(:,i) - ldata_all(end-N+1:end,i))./...
        abs(ldata_all(end-N+1:end,i)) ;
    figure(11);
    subplot(3,1,i);
    plot(1:N,err_rel(:,i))
    title(['Error relativo del ' nombre(i) 'en %'])
end

%% Bondad de ajuste de la predicción. Predicción estática.

RMSE = zeros(1,length(estimacion)); % Root Mean Squared Error
for i = 1:length(estimacion)
    RMSE(i) = sqrt(mean((ldata_all(end-N+1:end,i) - ldata_f(:,i)).^2));
end
%%%%%%%%RMSE%%%%%%%
% EURO-US : 0.2647
% Repsol  : 1.0659    
% Crude   : 1.0013
%%%%%%%%%%%%%%%%%%%

MAE = zeros(1,length(estimacion)); % Mean Absolute Error 
for i = 1:length(estimacion)
    MAE(i) = mean(abs(ldata_all(end-N+1:end,i) - ldata_f(:,i)));
end
%%%%%%%%MAE%%%%%%%%
% EURO-US : 0.1958
% Repsol  : 0.7600  
% Crude   : 0.7355
%%%%%%%%%%%%%%%%%%%

MAPE = zeros(1,length(estimacion)); % Mean Absolute Percentage Error 
for i = 1:length(estimacion)
    MAPE(i) = 100*mean(abs((ldata_all(end-N+1:end,i) - ldata_f(:,i)) ./...
        ldata_all(end-N+1:end,i)));
end
%%%%%%%%MAPE%%%%%%%
% EURO-US : 0.6538
% Repsol  : 0.2581    
% Crude   : 0.1598
%%%%%%%%%%%%%%%%%%%

TIC = zeros(1,length(estimacion)); % Theil Inequality Coefficient 
for i = 1:length(estimacion)
    TIC(i) = RMSE(i) / (mean(ldata_f(:,i).^2) * ...
        mean(ldata_all(end-N+1:end,i).^2));
end
%%%%%%%%TIC%%%%%%%%
% EURO-US : 0.3147
% Repsol  : 0.0001    
% Crude   : 0.0000
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SEGUNDA PARTE
% Datos menos N observaciones
data_N = data(1:end-N,:);
eurous = eurous(1:end-N); 
repsol = repsol(1:end-N);
crude = crude(1:end-N); 

%% Manipulación de los datos
% Estudiamos su ACF y PACF
figure(11);
subplot(3,2,1)
autocorr(eurous,20);
title('ACF del log del EURO-US');
subplot(3,2,2)
parcorr(eurous,20);
title('PACF del log del EURO-US');

subplot(3,2,3);
autocorr(repsol,20);
title('ACF del log de Repsol');
subplot(3,2,4);
parcorr(repsol,20);
title('PACF del log de Repsol');

subplot(3,2,5);
autocorr(crude,20);
title('ACF del log del crudo');
subplot(3,2,6);
parcorr(crude,20);
title('PACF del log del crudo');
% Necesitamos tomar diferencias: ACF decrece lenta y linealmente

%% Tomamos diferencias
deurous = eurous(2:end)-eurous(1:end-1);
drepsol = repsol(2:end)-repsol(1:end-1);
dcrude = crude(2:end)-crude(1:end-1);
ddata = data_N(2:end,:)-data_N(1:end-1,:);

figure(12);
subplot(3,2,1)
autocorr(deurous,20);
title('ACF del rendimiento del EURO-US');
subplot(3,2,2)
parcorr(deurous,20);
title('PACF del rendimiento del EURO-US');

subplot(3,2,3);
autocorr(drepsol,20);
title('ACF del rendimiento de Repsol');
subplot(3,2,4);
parcorr(drepsol,20);
title('PACF del rendimiento de Repsol');

subplot(3,2,5);
autocorr(dcrude,20);
title('ACF del rendimiento del crudo');
subplot(3,2,6);
parcorr(dcrude,20);
title('PACF del rendimiento del crudo');
%% Modelo
% Planteamos inicialmente un modelo ARIMA(1,1,1)

for i = 1:length(data(1,:))
    estimacion2(i) = estimate(arima(1,1,1),data_N(:,i));
end
% EURO-US y Repsol ningún parámetro significativo, Crude cte. no
% significativa

%% Análisis de los residuos

for i = 1:length(data(1,:))
    [res(:,i),varres(:,i),logL(i)] = infer(estimacion2(i),data_N(:,i));
    resstd(:,i) = res(:,i)/sqrt(estimacion2(i).Variance);
    [h(i),p(i),jbstat(i),critval(i)] = jbtest(resstd(:,i),alpha);
end
h, p, jbstat, critval
% Rechazamos la hipótesis de normalidad para todas las series

%% Gráficos de los residuos
nombre = ["EURO-US","Repsol","crudo"];
for i = 1:length(estimacion2)
    figure(12+i);
    subplot(2,2,1);
    plot(resstd(:,i));
    title(['Residuos estandarizados de', nombre(i)]);
    
    subplot(2,2,2);
    histogram(resstd(:,i),20);
    title(['Residuos estandarizados de ', nombre(i)]);
    
    subplot(2,2,3);
    autocorr(resstd(:,i),20);
    
    subplot(2,2,4);
    parcorr(resstd(:,i),20);
end

% Comprobación de los residuos

for i = 1:length(estimacion2)
    significativos_ma(:,i) = abs(autocorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ma(:,i) = IND'.*significativos_ma(2:end,i);
    significativos_ar(:,i) = abs(parcorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ar(:,i) = IND'.*significativos_ar(2:end,i);
end
lags_significativos_ar
lags_significativos_ma
% EURO-US: AR(6) y MA(6)
% Repsol : AR([3,5,17,18]) y MA([3,5,17,18])
% Crude  : AR([4,5,12,20]) y MA([4,5,12,20])

%% Reajustamos el modelo
% Eliminamos la constante en todos porque sale no significativa.
% EURO-US
estimacion2(1) = estimate(arima('Constant',0,'ARlags',6,'MAlags',6,...
                'D',1),data_N(:,1));
            
% Repsol
estimacion2(2) = estimate(arima('Constant',0,'ARlags',[3,5,17,18],'MAlags',...
                [3,5,17,18],'D',1),data_N(:,2));

% Crude            
estimacion2(3) = estimate(arima('Constant',0,'ARlags',[4,5,12,20],...
                'MAlags',[4,5,12,20],'D',1),data_N(:,3)); 

% EURO-US y Crude todos parámetros significativos. Repsol no significativos
% AR([3,18]) y MA([3,5,18]).
%% Volvemos a estudiar los residuos

for i = 1:length(data(1,:))
    [res(:,i),varres(:,i),logL(i)] = infer(estimacion2(i),data_N(:,i));
    resstd(:,i) = res(:,i)/sqrt(estimacion2(i).Variance);
    [h(i),p(i),jbstat(i),critval(i)] = jbtest(resstd(:,i),alpha);
end
h, p, jbstat, critval
% Rechazamos la hipótesis de normalidad para todas las series

%% Gráficas
for i = 1:length(estimacion2)
    figure(15+i);
    subplot(2,2,1);
    plot(resstd(:,i));
    title(['Residuos estandarizados de', nombre(i)]);
    
    subplot(2,2,2);
    histogram(resstd(:,i),20);
    title(['Residuos estandarizados de ', nombre(i)]);
    
    subplot(2,2,3);
    autocorr(resstd(:,i),20);
    
    subplot(2,2,4);
    parcorr(resstd(:,i),20);
end

% Comprobación de los residuos
for i = 1:length(estimacion2)
    significativos_ma(:,i) = abs(autocorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ma(:,i) = IND'.*significativos_ma(2:end,i);
    significativos_ar(:,i) = abs(parcorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ar(:,i) = IND'.*significativos_ar(2:end,i);
end
lags_significativos_ar
lags_significativos_ma

% EURO-US ruido blanco. Repsol significativo en AR([1,2]) y MA(1). 
% Crude en AR([1,2]) y MA([1,2]).

%% Reajustamos el modelo

% Repsol
estimacion2(2) = estimate(arima('Constant',0,'ARlags',[1,2,5,17],'MAlags',...
                [1,17],'D',1),data(:,2));

% Crude            
estimacion2(3) = estimate(arima('Constant',0,'ARlags',[1,2,4,5,12,20],...
                'MAlags',[1,2,4,5,12,20],'D',1),data(:,3)); 

for i = 1:length(estimacion2)
    significativos_ma(:,i) = abs(autocorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ma(:,i) = IND'.*significativos_ma(2:end,i);
    significativos_ar(:,i) = abs(parcorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ar(:,i) = IND'.*significativos_ar(2:end,i);
end
lags_significativos_ar
lags_significativos_ma

% Parámetros no significativos: 
% Repsol: MA(1)
% Crude : AR([1,2]) y MA([1,2]). 
% Respecto a los residuos, Repsol sigue siendo significativo en AR([1,2])
% y MA(1). Crudo en AR([1,2]) y MA([1,2]) a pesar de incluir los retardos.

%% Reajustamos el modelo

% Repsol
estimacion2(2) = estimate(arima('Constant',0,'ARlags',[1,2,5,17],'MAlags',...
                [17],'D',1),data(:,2));

% Crude            
estimacion2(3) = estimate(arima('Constant',0,'ARlags',[4,5,12,20],...
                'MAlags',[4,5,12,20],'D',1),data(:,3)); 

for i = 1:length(estimacion2)
    significativos_ma(:,i) = abs(autocorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ma(:,i) = IND'.*significativos_ma(2:end,i);
    significativos_ar(:,i) = abs(parcorr(resstd(:,i))) > banda_confianza;
    lags_significativos_ar(:,i) = IND'.*significativos_ar(2:end,i);
end
lags_significativos_ar
lags_significativos_ma

% Todos los parámetros son significativos excepto Repsol AR(2).
% Respecto a los residuos, Repsol sigue siendo significativo en AR([1,2])
% y MA(1). Crudo en AR([1,2]) y MA([1,2]) a pesar de incluir los retardos.

%% Último ajuste
% Repsol
estimacion2(2) = estimate(arima('Constant',0,'ARlags',[1,5,17],'MAlags',...
                [17],'D',1),data(:,2));
            
%% Predicción y bandas de confianza
% Predicción estática (a 1 día)
for i = 1:length(estimacion2)
    for j = 1:N
        [data_f(j,i),dt_errorf(j,i)] = forecast(estimacion2(i),1,'Y0',...
            data_N(:,i));
        b_inferior(j,i) = data_f(j,i) - Z*sqrt(dt_errorf(j,i));
        b_superior(j,i) = data_f(j,i) + Z*sqrt(dt_errorf(j,i));
        data_N = data(1:end-N+j,:);
    end
    data_N = data(1:end-N,:);
end
% Gráfica
M = round(1.5*N); % Puntos a mostrar
for i = 1:length(estimacion2)
    figure(19);
    subplot(3,1,i);
    plot(data_N(end-M:end,i),'Color',[.7,.7,.7]);
    hold on
    h1 = plot(length(data_N(end-M:end,i)):(N-1+length(data_N(end-M:end,i))),...
        b_inferior(:,i),'r:','LineWidth',2);
    plot(length(data_N(end-M:end,i)):(N-1+length(data_N(end-M:end,i))),...
        b_superior(:,i),'r:','LineWidth',2)
    h2 = plot(length(data_N(end-M:end,i)):(N-1+length(data_N(end-M:end,i))),...
        data_f(:,i),'k','LineWidth',2);
      legend([h1 h2], a, 'Predicción a 1 día', 'Location','NorthWest')
    title(['Predicción estática del ' nombre(i)])
    hold off
end

%% Predicción dinámica (a k días, con k = 1,2,...,N)
data_f2 = zeros(N,length(estimacion2)); 
dt_errorf2 = zeros(N,length(estimacion2));
for i = 1:length(estimacion2)
    [data_f2(:,i),dt_errorf2(:,i)] = forecast(estimacion2(i),N,'Y0',...
        data_N(:,i));
    b_inferior(:,i) = data_f2(:,i) - Z*sqrt(dt_errorf2(:,i));
    b_superior(:,i) = data_f2(:,i) + Z*sqrt(dt_errorf2(:,i));
end

for i = 1:length(estimacion2)
    figure(20);
    subplot(3,1,i);
    plot(data_N(end-M:end,i),'Color',[.7,.7,.7]);
    hold on
    h1 = plot(length(data_N(end-M:end,i)):(N-1+length(data_N(end-M:end,i))),...
        b_inferior(:,i),'r:','LineWidth',2);
    plot(length(data_N(end-M:end,i)):(N-1+length(data_N(end-M:end,i))),...
        b_superior(:,i),'r:','LineWidth',2)
    h2 = plot(length(data_N(end-M:end,i)):(N-1+length(data_N(end-M:end,i))),...
        data_f2(:,i),'k','LineWidth',2);
      legend([h1 h2], a, 'Predicción puntual', 'Location','NorthWest')
    title(['Predicción dinámica del ' nombre(i)])
    hold off
end

%% Comparación valores reales
for i = 1:length(estimacion2)
    figure(21);
    subplot(3,1,i);
    plot(data(end-M-N:end,i));
    title(['Valor real del ' nombre(i)])
end

% Error predicción
err_rel = zeros(N,length(estimacion2));
for i = 1:length(estimacion2)
    err_rel(:,i) = 100*abs(data_f(:,i) - data(end-N+1:end,i))./...
        abs(data(end-N+1:end,i)) ;
    figure(11);
    subplot(3,1,i);
    plot(1:N,err_rel(:,i))
    title(['Error relativo del ' nombre(i) 'en %'])
end

%% Bondad de ajuste de la predicción. Predicción estática.

RMSE2 = zeros(1,length(estimacion2)); % Root Mean Squared Error
for i = 1:length(estimacion2)
    RMSE2(i) = sqrt(mean((data(end-N+1:end,i) - data_f(:,i)).^2));
end
%%%%%%%%RMSE%%%%%%%
% EURO-US : 0.0020
% Repsol  : 0.2006    
% Crude   : 0.9789
%%%%%%%%%%%%%%%%%%%

MAE2 = zeros(1,length(estimacion2)); % Mean Absolute Error 
for i = 1:length(estimacion2)
    MAE2(i) = mean(abs(data(end-N+1:end,i) - data_f(:,i)));
end
%%%%%%%%MAE%%%%%%%%
% EURO-US : 0.0014
% Repsol  : 0.1417    
% Crude   : 0.7202
%%%%%%%%%%%%%%%%%%%

MAPE2 = zeros(1,length(estimacion)); % Mean Absolute Percentage Error 
for i = 1:length(estimacion)
    MAPE2(i) = 100*mean(abs((data(end-N+1:end,i) - data_f(:,i)) ./...
        data(end-N+1:end,i)));
end
%%%%%%%%MAPE%%%%%%%
% EURO-US : 0.1955
% Repsol  : 0.7453    
% Crude   : 0.7212
%%%%%%%%%%%%%%%%%%%

TIC2 = zeros(1,length(estimacion)); % Theil Inequality Coefficient 
for i = 1:length(estimacion)
    TIC2(i) = RMSE(i) / (mean(data_f(:,i).^2) * ...
        mean(data(end-N+1:end,i).^2));
end
%%%%%%%%TIC%%%%%%%
% EURO-US : 0.0066
% Repsol  : 0.0000    
% Crude   : 0.0000
%%%%%%%%%%%%%%%%%%%

%% Diebold y Mariano

error1 = ldata_all(end-N+1:end,i) - ldata_f(:,i);
error2 = data(end-N+1:end,i) - data_f(:,i);
d = abs(error1)-abs(error2);
dbar = mean(d);
gamma0 = var(d);
N_cov = round(N^(1/3)-0.5)+1;
[acf,lags,bounds] = autocorr(d,N_cov);
vgamman = acf(2:end-1,1) * gamma0;
var_dbar = (1/N)*(gamma0 + 2*sum(vgamman));

DM = dbar/(var_dbar^0.5)
p_value = 2*(1-normcdf(abs(DM),0,1))
% Rechazamos la H_0: esto quiere decir que los modelos no tienen la misma
% capacidad predictiva.
