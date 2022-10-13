%+------------------------------------------------------------------------+
%    Ricardo Vega, Vladimir Rabinovitch Likhtman, Raúl Castillo Pérez
%                     Método SPPS para Ecuación de Dispersión
%+------------------------------------------------------------------------+

% clc;
clear all;
close all;

%% Raíces de la ecuación de dispersión

H = 182.09; % Alturas (Las cuales se puden modificar para constatar los NMO obtenidos entre
               % el método SPPS y el  analítico)
w = 83E6; % Frecuencia con la que se operará
c0 = 3E8;
wp = 2*pi*4.429E6; %Frecuencias del plasma en respectivas H
ka = w/c0; % Numero de onda en la atmósfera

ki = sqrt(w^2-wp^2)/c0; % Numero de onda en la ionósfera
alfa = linspace(ki,ka,100);
alfa = alfa(2:end-1);
eqn = tan(sqrt(ka^2 - alfa.^2)*H) + sqrt(ka^2 - alfa.^2)./sqrt(alfa.^2 - ki^2);


%% Graficación del método analítico
Eqn = spline(alfa, eqn); % Construye una función continua aproximada con valores discretos dados
NMO = fnzeros(Eqn); %localiza los ceros del spline
NMO = NMO(1,:)
figure
plot(alfa, eqn)
grid on
title('\Lambda(\omega,\alpha)')
xlabel('\alpha')
legend('\Lambda')