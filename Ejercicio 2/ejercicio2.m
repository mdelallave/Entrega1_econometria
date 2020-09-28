%% Limpiamos
close('all')
clearvars
clc

%% Importar datos
data = xlsread('data.csv','B4:D6613');
eurous = data(:,1); % Tipo de cambio euro-us
repsol = data(:,2); % Repsol
crude = data(:,3);  % Futuros del crudo

%% Manipulación de los datos
% Transformación logarítmica (rendimientos)
leurous = 100*log(eurous); lrepsol = 100*log(repsol); 
lcrude = 100*log(crude); ldata = 100*log(data);
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

% Tomamos diferencias
dleurous = leurous(2:end)-leurous(1:end-1);
dlrepsol = lrepsol(2:end)-lrepsol(1:end-1);
dlcrude = lcrude(2:end)-lcrude(1:end-1);
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
pause
%% Modelo

for i = 1:length(data(1,:))
    estimacion(i) = estimate(arima(1,1,1),ldata(:,i));
end
pause
%%
% Residuos y test de Jarque-Bera
res = zeros(length(data(:,1)),length(data(1,:)));
varres = zeros(length(data(:,1)),length(data(1,:)));
logL = zeros(length(data(1,:)),1);
resstd = zeros(length(data(:,1)),length(data(1,:)));
h = zeros(length(data(1,:)),1);
p = zeros(length(data(1,:)),1);
jbstat = zeros(length(data(1,:)),1);
critval = zeros(length(data(1,:)),1);

for i = 1:length(data(1,:))
    [res(:,i),varres(:,i),logL(i)] = infer(estimacion(i),ldata(:,i));
    resstd(:,i) = res(:,i)/sqrt(estimacion(i).Variance);
    [h(i),p(i),jbstat(i),critval(i)] = jbtest(resstd(:,i),0.01);
end
h, p, jbstat, critval
