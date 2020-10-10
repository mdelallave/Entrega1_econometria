%% Limpiamos
close('all')
clearvars
clc

%% Importar datos
data = xlsread('data.csv','B4:D6613');
eurous = data(:,1); % Tipo de cambio euro-us
repsol = data(:,2); % Repsol
crude = data(:,3);  % Futuros del crudo

%% Par�metros iniciales
N = 100; % Cantidad de datos omitidos.
alpha = 0.05; % Nivel de significaci�n.
T = length(data(1:end-N,1)); % Cantidad de datos.
%% Manipulaci�n de los datos
% Transformaci�n logar�tmica (rendimientos)
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

%% Tomamos diferencias
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
%% Selecci�n modelo
pp = 5; qq = 5; % Maxima p y q que dejamos que el modelo tome.
D = 1; % N�mero de diferencias de nuestro modelo

% Tipo de cambio EURO/US
check_eurous = checkArima(leurous,pp,D,qq,T);
[M,I] = min(check_eurous(:)); % I: �ndice que contiene el valor m�s peque�o
[p, q] = ind2sub(size(check_eurous),I); % Coordenadas del valor m�s peque�o

model_eurous = estimate(arima(p,D,q),leurous);
pause

% Repsol
check_repsol = checkArima(lrepsol,pp,D,qq,T);
[M,I] = min(check_repsol(:)); 
[p, q] = ind2sub(size(check_repsol),I);

model_repsol = estimate(arima(p,D,q),lrepsol);
pause

% Crudo
check_crude = checkArima(lcrude,pp,D,qq,T);
[M,I] = min(check_crude(:)); 
[p, q] = ind2sub(size(check_crude),I); 

model_crude = estimate(arima(p,D,q),lcrude);
pause
%%
for i = 1:length(ldata(1,:))
    estimacion(i) = estimate(arima(1,1,1),ldata(:,i));
end
%%
% Residuos y test de Jarque-Bera
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