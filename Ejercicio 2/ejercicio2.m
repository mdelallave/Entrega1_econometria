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
nombre = ["EURO-US","Repsol","crudo"]; % Títulos de los gráficos

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

for i = 1:length(dldata(1,:))
   figure(3);
   subplot(3,1,i);
   plot(dldata(:,i))
   title(['Diferencia log de ', nombre(i)])
end
% Podemos apreciar que nuestras series son estacionales y no hace falta que
% tomemos segundas diferencias.

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
for i = 1:length(estimacion)
    figure(3+i);
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
    figure(6+i);
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
    figure(10);
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
    figure(11);
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
    figure(12);
    subplot(3,1,i);
    plot(ldata_all(end-M-N:end,i));
    title(['Valor real del ' nombre(i)])
end

% Error predicción
err_rel = zeros(N,length(estimacion));
for i = 1:length(estimacion)
    err_rel(:,i) = 100*abs(ldata_f(:,i) - ldata_all(end-N+1:end,i))./...
        abs(ldata_all(end-N+1:end,i)) ;
    figure(13);
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
dataN = data(1:end-N,:); eurousN = eurous(1:end-N); 

% Número de lags de cada modelo
for i = 1:length(estimacion)
    P(i) = length(estimacion(i).AR);
    Q(i) = length(estimacion(i).MA);
end
max_lags = max(max(P,Q));

% Inicialización de variables
const = zeros(1,length(estimacion));
ar_parameters = zeros(max(P),length(estimacion));
ma_parameters = zeros(max(Q),length(estimacion));
var = zeros(1,length(estimacion));
epsilon = zeros(N,length(estimacion));

% Expresión matricial de todos los parámetros de nuestros modelos y de la
% varianza. Cada columna corresponde a una serie temporal.
for i = 1:length(estimacion)
    const(i) = estimacion(i).Constant;
    ar = estimacion(i).AR;
    ma = estimacion(i).MA;
    ar = cell2mat(ar);
    ma = cell2mat(ma);
    ar_parameters(:,i) = [ar(1:end), zeros(1,max(P)-length(ar))];
    ma_parameters(:,i) = [ma(1:end), zeros(1,max(Q)-length(ma))];
    var(i) = estimacion(i).Variance;
    epsilon(:,i) = normrnd(0,sqrt(var(i)),1,N);
    [res(:,i),varres(:,i),logL(i)] = infer(estimacion(i),ldata(:,i));
end
ar_parameters = flip(ar_parameters)/100; % Le damos la vuelta
ma_parameters = flip(ma_parameters)/100;

% A los residuos que ya teníamos, le añadimos las siguientes 100 
% innovaciones N~(0,var(i))
eps_total = [res; epsilon]; 

%% Predicción y bandas de confianza
%Inicialización de variables
data_f = zeros(N + max_lags + 1, length(estimacion));
b_inferior_precios = zeros(N,length(estimacion));
b_superior_precios = zeros(N,length(estimacion));

for i = 1:length(estimacion)
    data_f(1:max_lags+1,i) = dataN(end-max_lags:end,i);
    for j = max_lags+1:N+max_lags
        data_f(j+1,i) = exp((log(data_f(j,i)) + 1/100 * (const(i) + ...
            ar_parameters(:,i)' * (log(data_f(j-max_lags+1:j,i)) - ...
            log(data_f(j-max_lags:j-1,i))) - ma_parameters(:,i)' * ...
            eps_total(end-N-2*max_lags+j:end-N-max_lags+j-1,i)))+1/2e4*var(i));
        b_inferior_precios(j-max_lags,i) = data_f(j+1,i) - Z*sqrt(var(i));
        b_superior_precios(j-max_lags,i) = data_f(j+1,i) + Z*sqrt(var(i));
    end
end
% Gráfica
M = round(1.5*N); % Puntos a mostrar
for i = 1:length(estimacion)
    figure(14);
    subplot(3,1,i);
    plot(dataN(end-M:end,i),'Color',[.7,.7,.7]);
    hold on
    h1 = plot(length(dataN(end-M:end,i)):(N-1+length(dataN(end-M:end,i))),...
        b_inferior_precios(:,i),'r:','LineWidth',2);
    plot(length(dataN(end-M:end,i)):(N-1+length(dataN(end-M:end,i))),...
        b_superior_precios(:,i),'r:','LineWidth',2)
    h2 = plot(length(dataN(end-M:end,i)):(N-1+length(dataN(end-M:end,i))),...
        data_f(max_lags+2:end,i),'k','LineWidth',2);
      legend([h1 h2], a, 'Predicción a 1 día', 'Location','NorthWest')
    title(['Predicción estática del ' nombre(i)])
    hold off
end

%% Comparación valores reales
for i = 1:length(estimacion)
    figure(15);
    subplot(3,1,i);
    plot(data(end-M-N:end,i));
    title(['Valor real del ' nombre(i)])
end

%% Error predicción
err_rel = zeros(N,length(estimacion));
for i = 1:length(estimacion)
    err_rel(:,i) = 100*abs(data_f(max_lags+2:end,i) - ...
        data(end-N+1:end,i))./ abs(data(end-N+1:end,i)) ;
    figure(16);
    subplot(3,1,i);
    plot(1:N,err_rel(:,i))
    title(['Error relativo del ' nombre(i) 'en %'])
end

%% Bondad de ajuste de la predicción (precios). Predicción estática.

RMSE2 = zeros(1,length(estimacion)); % Root Mean Squared Error
for i = 1:length(estimacion)
    RMSE2(i) = sqrt(mean((data(end-N+1:end,i) - ...
        data_f(max_lags+2:end,i)).^2));
end
%%%%%%%%RMSE%%%%%%%
% EURO-US : 0.0182
% Repsol  : 0.7094   
% Crude   : 6.2105
%%%%%%%%%%%%%%%%%%%

MAE2 = zeros(1,length(estimacion)); % Mean Absolute Error 
for i = 1:length(estimacion)
    MAE2(i) = mean(abs(data(end-N+1:end,i) - data_f(max_lags+2:end,i)));
end
%%%%%%%%MAE%%%%%%%%
% EURO-US : 0.0146
% Repsol  : 0.5639   
% Crude   : 4.4280
%%%%%%%%%%%%%%%%%%%

MAPE2 = zeros(1,length(estimacion)); % Mean Absolute Percentage Error 
for i = 1:length(estimacion)
    MAPE2(i) = 100*mean(abs((data(end-N+1:end,i) - ...
        data_f(max_lags+2:end,i)) ./ data(end-N+1:end,i)));
end
%%%%%%%%MAPE%%%%%%%
% EURO-US : 1.9637
% Repsol  : 0.7453    
% Crude   : 0.7212
%%%%%%%%%%%%%%%%%%%

TIC2 = zeros(1,length(estimacion)); % Theil Inequality Coefficient 
for i = 1:length(estimacion)
    TIC2(i) = RMSE(i) / (mean(data_f(max_lags+2:end,i).^2) * ...
        mean(data(end-N+1:end,i).^2));
end
%%%%%%%%TIC%%%%%%%%
% EURO-US : 0.9193
% Repsol  : 8 e-06  
% Crude   : 8 e-09
%%%%%%%%%%%%%%%%%%%

%% Vamos a comparar los estadísticos:

RMSE_comparacion = RMSE > RMSE2
MAE_comparacion = MAE > MAE2
MAPE_comparacion = MAPE > MAPE2
TIC_comparacion = TIC > TIC2

% Tanto el RMSE como el MAE de los rendimientos es mayor que si tratamos
% las series en precios excepto para el Crude en el MAE. Sin embargo es al
% revés si usamos el MAPE y el TIC.

% Dado que el RMSE mide la variabilidad que tienen los errores de previsión,
% podemos decir que la predicción del modelo cuando usamos precios es más 
% acertada que cuando usamos rendimientos.
