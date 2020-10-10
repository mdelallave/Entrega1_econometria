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
    [h(i),p(i),jbstat(i),critval(i)] = jbtest(resstd(:,i),0.01);
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
% Residuos EURO-US son ruido blanco, pero para Repsol y el crudo hay que
% ajustar el modelo. Todos los coeficientes de EURO-US son no
% significativos, por lo que hay que especificar otro tipo de modelo.

%% Reajustamos el modelo

% EURO-US
estimacion(1) = estimate(arima('Constant',0,'ARlags',6,'MAlags',1,...
                'D',1),ldata(:,1));
            
% Repsol
estimacion(2) = estimate(arima('Constant',0,'ARlags',1,'MAlags',[3,5],...
                'D',1),ldata(:,2));

% Crude
estimacion(3) = estimate(arima('Constant',0,'ARlags',[2,4],'MAlags',1,...
                'D',1),ldata(:,3));
            
% estimacion(3) = estimate(arima('Constant',0,'ARlags',[2,4],'MAlags',[2,4],...
%                 'D',1),ldata(:,3));
            % Así sí que sale todo significativo.
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
    [h(i),p(i),jbstat(i),critval(i)] = jbtest(resstd(:,i),0.01);
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
% Los residuos de EUROUS ahora se comportan como un r.b. Repsol falla en el
% lag 2 y el crudo en el 3 y el 5.


%% Predicción y bandas de confianza
% Inicialización de variables
ldata_f = zeros(N,length(estimacion)); 
dt_errorf = zeros(N,length(estimacion));
b_inferior = zeros(N,length(estimacion));
b_superior = zeros(N,length(estimacion));

for i = 1:length(estimacion)
    [ldata_f(:,i),dt_errorf(:,i)] = forecast(estimacion(i),N,'Y0',ldata(:,i));
    b_inferior(:,i) = ldata_f(:,i) - norminv(alpha/2,0,1)*sqrt(dt_errorf(:,i));
    b_superior(:,i) = ldata_f(:,i) + norminv(alpha/2,0,1)*sqrt(dt_errorf(:,i));
end
% Gráficas
for i = 1:length(estimacion)
    figure(8);
    subplot(3,1,i);
    plot(ldata(:,i),'Color',[.7,.7,.7]);
    hold on
    h1 = plot(length(ldata(:,i)):(N-1+length(ldata(:,i))),b_inferior(:,i),...
        'r:','LineWidth',2);
    plot(length(ldata(:,i)):(N-1+length(ldata(:,i))),b_superior(:,i),...
        'r:','LineWidth',2)
    h2 = plot(length(ldata(:,i)):(N-1+length(ldata(:,i))),ldata_f(:,i),...
        'k','LineWidth',2);
      legend([h1 h2], a, 'Predicción puntual', 'Location','NorthWest')
    title(['Predicción del ' nombre(i)])
    hold off
end
