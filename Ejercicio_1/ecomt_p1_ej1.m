%% Limpiamos
clear all
close all
clc

%% Muestra
M = 1000;
T = [20,50,100];
mu = 0;
sigma = 1;
signif_level = [0.01, 0.05, 0.1];
X1 = normrnd(mu,sigma,M,T(1));
X2 = normrnd(mu,sigma,M,T(2));
X3 = normrnd(mu,sigma,M,T(3));

JB_X1= zeros(M,length(signif_level));
JB_X2= zeros(M,length(signif_level));
JB_X3= zeros(M,length(signif_level));

for j = 1:length(signif_level)
    for i = 1:M
        JB_X1(i,j) = jbtest(X1(i,:),signif_level(j));
        JB_X2(i,j) = jbtest(X2(i,:),signif_level(j));
        JB_X3(i,j) = jbtest(X3(i,:),signif_level(j));
    end
end

n_rejects_X1 = sum(JB_X1 == 1);
per_rejects_X1 = n_rejects_X1/M*100;
n_rejects_X2 = sum(JB_X2 == 1);
per_rejects_X2 = n_rejects_X2/M*100;
n_rejects_X3 = sum(JB_X3 == 1);
per_rejects_X3 = n_rejects_X3/M*100;


%%%%% Probamos a calcular por nosotros mismos el estadístico Bera-Jarque
% sk_X1 = skewness(X1(1,:),1); % No corregir el sesgo
% kur_X1 = kurtosis(X1(1,:),1); % No corregir el sesgo
% 
% myJB_X1 = T(1)/6*(sk_X1.^2 + 0.25*(kur_X1-3).^2) % Nuestro calculo
% [h,p,jbstat,critval] = jbtest(X1(1,:),0.05) % Matlab's function
% % Results:
% % myJB_X1 =
% % 
% %     1.4387
% % jbstat =
% % 
% %     1.4387
%%%%%

%% Contraste de potencia
Z1 = zeros(M,T(1),30);
Z2 = zeros(M,T(2),30);
Z3 = zeros(M,T(3),30);
JB_Z1= zeros(M,length(signif_level),30);
JB_Z2= zeros(M,length(signif_level),30);
JB_Z3= zeros(M,length(signif_level),30);
n_rejects_Z1 = zeros(30,3);
per_rejects_Z1 = zeros(30,3);
n_rejects_Z2 = zeros(30,3);
per_rejects_Z2 = zeros(30,3);
n_rejects_Z3 = zeros(30,3);
per_rejects_Z3 = zeros(30,3);

for n = 1:30
    Z1(:,:,n) = trnd(n,M,T(1));
    Z2(:,:,n) = trnd(n,M,T(2));
    Z3(:,:,n) = trnd(n,M,T(3));
    
    for j = 1:length(signif_level)
        for i = 1:M
            JB_Z1(i,j,n) = jbtest(Z1(i,:,n),signif_level(j));
            JB_Z2(i,j,n) = jbtest(Z2(i,:,n),signif_level(j));
            JB_Z3(i,j,n) = jbtest(Z3(i,:,n),signif_level(j));
        end
    end
    
    n_rejects_Z1(n,:) = sum(JB_Z1(:,:,n) == 1,1);
    per_rejects_Z1(n,:) = n_rejects_Z1(n,:)/M*100;
    n_rejects_Z2(n,:) = sum(JB_Z2(:,:,n) == 1,1);
    per_rejects_Z2(n,:) = n_rejects_Z2(n,:)/M*100;
    n_rejects_Z3(n,:) = sum(JB_Z3(:,:,n) == 1,1);
    per_rejects_Z3(n,:) = n_rejects_Z3(n,:)/M*100;
end

figure(1);
plot(1:30,per_rejects_Z1(:,1),'.-',1:30,per_rejects_Z1(:,2),'-o',...
          1:30,per_rejects_Z1(:,3),'-*')
%title('Potencia del contraste Bera-Jarque (T = 20)')
legend('$\alpha = 1\%$','$\alpha = 5\%$','$\alpha = 10\%$','Interpreter','latex')
xlabel('Grados de libertad t-Student')
ylabel('Porcentaje de rechazos de la Ho')

figure(2);
plot(1:30,per_rejects_Z2(:,1),'.-',1:30,per_rejects_Z2(:,2),'-o',...
          1:30,per_rejects_Z2(:,3),'-*')
%title('Potencia del contraste Bera-Jarque (T = 50)')
legend('$\alpha = 1\%$','$\alpha = 5\%$','$\alpha = 10\%$','Interpreter','latex')
xlabel('Grados de libertad t-Student')
ylabel('Porcentaje de rechazos de la Ho')

figure(3);
plot(1:30,per_rejects_Z3(:,1),'.-',1:30,per_rejects_Z3(:,2),'-o',...
          1:30,per_rejects_Z3(:,3),'-*')
%title('Potencia del contraste Bera-Jarque (T = 100)')
legend('$\alpha = 1\%$','$\alpha = 5\%$','$\alpha = 10\%$','Interpreter','latex')
xlabel('Grados de libertad t-Student')
ylabel('Porcentaje de rechazos de la Ho')
