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
% Estudiamos su estacionariedad
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
% No estacionario, necesitamos tomar diferencias.

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
significativos = zeros(21,length(estimacion));
lags_significativos = zeros(20,length(estimacion));
IND = 1:20; % Índice

for i = 1:length(estimacion)
    significativos(:,i) = abs(autocorr(resstd(:,i))) > banda_confianza;
    lags_significativos(:,i) = IND'.*significativos(2:end,i);
end
lags_significativos

% EURO-US: 6
% Repsol : [2,3,5]
% Crude  : [1,2,4,8,15]

%% Reajustamos el modelo

% EURO-US
estimacion(1) = estimate(arima('Constant',0,'ARlags',6,'MAlags',6,...
                'D',1),ldata(:,1));
            
% Repsol
estimacion(2) = estimate(arima('Constant',0,'ARlags',[2,3,5],'MAlags',1,...
                'D',1),ldata(:,2));

% Crude            
estimacion(3) = estimate(arima('Constant',0,'ARlags',[1,2,4,8,15],'MAlags',...
                1,'D',1),ldata(:,3)); 
            
% Así sí que sale todo significativo. Si solo ponemos los lags en AR,
% entonces MA(1) de EURO-Us  no salen significativos. Si ponemos los lags
% en AR y MA, la estimación del crudo (3) sale fatal. En Repsol sale bien
% de las dos formas.
                       
%% Volvemos a estudiar los residuos
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
    significativos(:,i) = abs(autocorr(resstd(:,i))) > banda_confianza;
    lags_significativos(:,i) = IND'.*significativos(2:end,i);
end
lags_significativos

% EURO-US y Repsol ruido blanco. Crudo significativo en 5 y 6.

%% Reajustamos el modelo

% AR(4) de Crude no significativo            
estimacion(3) = estimate(arima('Constant',0,'ARlags',[1,2,5,6,8,15],'MAlags',...
                1,'D',1),ldata(:,3)); 

for i = 1:length(estimacion)
    significativos(:,i) = abs(autocorr(resstd(:,i))) > banda_confianza;
    lags_significativos(:,i) = IND'.*significativos(2:end,i);
end
lags_significativos

% Crudo sigue siendo significativo en 5 y 6. AIC menor y BIC mayor, por lo
% que volvemos al modelo anterior pero quitando AR(4) que no es
% significativo.

%%
estimacion(3) = estimate(arima('Constant',0,'ARlags',[1,2,8,15],'MAlags',...
                1,'D',1),ldata(:,3)); 

            % Mejoran tanto el AIC como el BIC eliminano AR(4)
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

MAE = zeros(1,length(estimacion)); % Mean Absolute Error 
for i = 1:length(estimacion)
    MAE(i) = mean(abs(ldata_all(end-N+1:end,i) - ldata_f(:,i)));
end

MAPE = zeros(1,length(estimacion)); % Mean Absolute Percentage Error 
for i = 1:length(estimacion)
    MAPE(i) = 100*mean(abs((ldata_all(end-N+1:end,i) - ldata_f(:,i)) ./...
        ldata_all(end-N+1:end,i)));
end

TIC = zeros(1,length(estimacion)); % Theil Inequality Coefficient 
for i = 1:length(estimacion)
    TIC(i) = RMSE(i) / (mean(ldata_f(:,i).^2) * ...
        mean(ldata_all(end-N+1:end,i).^2));
end
