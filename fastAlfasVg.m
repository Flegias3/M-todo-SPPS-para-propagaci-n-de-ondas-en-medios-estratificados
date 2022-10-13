%+------------------------------------------------------------------------+
%    Ricardo Vega, Vladimir Rabinovitch Likhtman, Raúl Castillo Pérez
%                     Método SPPS para Velocidades de grupo
%+------------------------------------------------------------------------+

clc;
close all;
clear all;
load awj.mat
format long

% Variables a considerar para obtener las soluciones
w = 83E6:13E3:95.987E6; %Frecuencia de operación
% awj = [];
% 
% for a = 1: length(w)
% c0 = 3E8; %Velocidad de la luz en el espacio libre.
% H = 81.03; %Alturas
% H1a = 105;
% H1b = 182.09;
% H1c = 255;
% H2 = 300;
% H1_H = H1c - H;%Intervalo de integración
% wpa = 2*pi*1.5457E6; %Frecuencia de plasma (ionósfera)
% wpb = 2*pi*3.148E6;
% wpc = 2*pi*4.429E6;
% wp1 = 2*pi*7E6; %Frecuencia de plasma (sobre la ionósfera)
% k1 = w(a)/c0; %Número de modo onda en atmósfera 
% k2a = sqrt(w(a)^2-wpa^2)/c0; %Número de modo onda en la ionósfera
% k2b = sqrt(w(a)^2-wpb^2)/c0;
% k2c = sqrt(w(a)^2-wpc^2)/c0;
% k3 = sqrt(w(a)^2-wp1^2)/c0; %Número de modo onda sobre la ionósfera
% dim = 1000; %Dimensión vectorial
% N = 50; %Número de potencias a considerar
% x0 = linspace(0,H,dim);
% x = linspace(H, H1c, dim);
% x3 = linspace(H1c,H2,dim);
% alfa = linspace(k3, k1, dim+2);
% alfa = alfa(2:end-1);
% 
% %------------------Para el caso 0 <= z <= H-----------------
% fi1 = sin(sqrt(k1^2-alfa.^2)'*x0);
% fi1Der = cos(sqrt(k1^2-alfa.^2)'*x0).*sqrt(k1^2-alfa.^2)';
% 
% %% Para el caso H <= z <= H1
% 
% % Solución de la ecuación homogenea
% p = ones(1, dim);
% q = zeros(1, dim);
% xDensity = dim/H1_H;
% r = -1*k2a^2.*ones(1, round(xDensity*(H1a-H)));
% r = [r -1*k2b^2.*ones(1, round(xDensity*(H1b-H1a)))];
% r = [r -1*k2c^2.*ones(1, round(xDensity*(H1c-H1b)))];
% lambda = 1;
% W = 1;
% WDer = 0;
% 
% % Construcción de potencias formales
% Xvirg = ones(1, dim);
% X = ones(1, dim);
% 
% % Para "v"
% for k = 1: N
%     if rem(k, 2) == 1
%         Xvirg(k+1, :) = ninteg(Xvirg(k, :)*W^2.*r, H1_H);
%         X(k+1, :) = ninteg(X(k, :)./(W^2.*p), H1_H);
%     else
%         X(k+1, :) = ninteg(X(k, :)*W^2.*r, H1_H);
%         Xvirg(k+1, :) = ninteg(Xvirg(k, :)./(W^2.*p), H1_H);
%     end
% end
% 
% % ----------Construcción de soluciones-------------------
% v1 = zeros(1, dim);
% v2 = zeros(1, dim);
% v1Der = zeros(1, dim);
% v2Der = zeros(1, dim);
% 
% for k = 0: N/2
%     v1 = v1 + W*lambda.*Xvirg(2*k+1, :);
%     if k > 0
%         v1Der = v1Der + lambda*Xvirg(2*k, :)./(W*p);
%     end
% end                                       
% 
% for k = 0: N/2-1
%     v2 = v2 + lambda*W.*X(2*k+2, :);
%     v2Der = v2Der + lambda*X(2*k+1, :)./(W*p);
% end
% 
% % Para "u"
% v = v1 + 1i*v2;
% vDer = v1Der + 1i*v2Der;
% 
% p = ones(1, dim);
% q = -r;
% %q = k2^2.*ones(1, dim);
% r = ones(1, dim);
% 
% for k = 1: N 
%     if rem(k, 2) == 1
%         Xvirg(k+1, :) = ninteg(Xvirg(k, :).*v.^2.*r, H1_H);
%         X(k+1, :) = ninteg(X(k, :)./(v.^2.*p), H1_H);
%     else
%         X(k+1, :) = ninteg(X(k, :).*v.^2.*r, H1_H);
%         Xvirg(k+1, :) = ninteg(Xvirg(k, :)./(v.^2.*p), H1_H);
%     end
% end
% 
% %--------------------------------------------------------------------------
% %Construcción de la solución en 50 km < = z <= 250 km
% u1 = zeros(N/2+1,dim);
% u2 = zeros(N/2+1,dim);
% u1Der = zeros(N/2+1,dim);
% u2Der = zeros(N/2+1,dim);
% 
% for cont = 0:N/2
%     u1(cont+1,:) = v.*Xvirg(1+2*cont,:);
%     if cont > 0
%         u1Der(cont+1,:) = Xvirg(2*cont,:)./(v.*p);
%     end
% end
% 
% u1Der = u1.*vDer./v + u1Der;
% u1 = flip(u1);
% u1Der = flip(u1Der);
% 
% for cont = 0:N/2-1
%     u2(cont+1,:) = v.*X(2+2*cont,:);
%     u2Der(cont+1,:) = X(1+2*cont,:)./(v.*p);
% end
% 
% u2Der = u2.*vDer./v + u2Der;
% u2 = flip(u2);
% u2Der = flip(u2Der);
% 
% %----------------Construcción de fi2 y fi2Der---------------------------
% u1Val = zeros(dim,dim);
% u2Val = zeros(dim,dim);
% u1DerVal = zeros(dim,dim);
% u2DerVal = zeros(dim,dim);
% 
% for cont = 1:dim        %Evaluación para cada punto en x (columnas)
%     u1Val(:,cont) = polyval(u1(:,cont),alfa.^2); %renglones con valores para alfa
%     u2Val(:,cont) = polyval(u2(:,cont),alfa.^2);
%     u1DerVal(:,cont) = polyval(u1Der(:,cont),alfa.^2);
%     u2DerVal(:,cont) = polyval(u2Der(:,cont),alfa.^2);
% end
% 
% C1 = fi1(:,end)./u1Val(:,1);
% C2 = fi1Der(:,end)./u2DerVal(:,1);
% 
% fi2 = real(u1Val).*C1 + real(u2Val).*C2;
% fi2Der = real(u1DerVal).*C1 + real(u2DerVal).*C2;
% 
% fi2A= fi2(:, end);
% fi2DerA = fi2Der(:, end);
% 
% %-------------------Ecuación de Dispersión-----------------------------
% EcDisp = fi2(:, end) + fi2Der(:, end)./sqrt(alfa.^2-k3^2)';
% EcDispSpline = spline(alfa.^2,EcDisp);
% Eigenvalores = fnzeros(EcDispSpline);
% NMO = sqrt(Eigenvalores(1,:));
% awj(a,:) = NMO;
% end

%-------------------Ecuación de Dispersión para derivar-----------------------------
% for b = 1: length(awj(1, :))
%     fi2Spline = spline(awj(:, b).^2, fi2(:, end));
%     fi2DerSpline = spline(awj(:, b).^2, fi2Der(:, end));
%     fi2awj(:, b) = fnval(fi2Spline, awj(:, b));
%     fi2Derawj(:, b) = fnval(fi2DerSpline, awj(:, b));
% end
% 
% for c = 1: length(awj(1, :))
%     fi2wSpline = spline(w, fi2awj(:, c));
%     fi2DerwSpline = spline(w, fi2Derawj(:, c));
%     fi2w(:, c) = fnval(w, fnder(fi2wSpline));
%     fi2Derw(:, c) = fnval(w, fnder(fi2DerwSpline));
%     fi2aSpline = spline(awj(:, c), fi2awj(:, c));
%     fi2DeraSpline = spline(awj(:, c), fi2Derawj(:, c));
%     fi2a(:, c) = fnval(awj(:, c), fnder(fi2aSpline));
%     fi2Dera(:, c) = fnval(awj(:, c), fnder(fi2DeraSpline));
% end
% 
% k3w = spline(w, k3);
% k3w = fnval(w, fnder(k3w));
% 
% for d = 1: length(awj(1,:))
%     Lw(:, d) = fi2w(:, d) - (k3'.*k3w'./sqrt(awj(:, d).^2 - k3'.^2)).*fi2Derawj(:, d) + ...
%         fi2Derw(:, d)./sqrt(awj(:, d).^2 - k3'.^2);
%     La(:, d) = fi2a(:, d) + (awj(:, d)./sqrt(awj(:, d).^2 - k3'.^2)).*fi2Derawj(:, d) + ...
%         fi2Dera(:, d)./sqrt(awj(:, d).^2 - k3'.^2);
% end
% 
% for e = 1: length(awj(1,:))
%     Vg(:, e) =  La(:, e)./Lw(:, e);
% end

%% -------------------Componentes para las Vg-----------------------------

% Graficas de las funciones de los numeros de modo de onda

for b = 1: length(awj(1,:))
    awSpline = spline(w, awj(:, b));
    awp(:, b) = fnval(fnder(awSpline), w);
end

for c =1: length(awj(1,:))
    Vg(:, c) = 1./awp(:, c);
end

figure
plot(w, awj(:,1), 'r')
hold on
grid on
plot(w, awj(:,2), 'g')
plot(w, awj(:,3), 'b')
title('Números de modo de onda (\alpha_{j}(\omega))')
xlabel('\omega (rad/s)')
ylabel('\alpha_{j}(\omega) (1/m)')
legend('\alpha_{1}(\omega)', '\alpha_{2}(\omega)', '\alpha_{3}(\omega)')
hold off

figure
plot(w, Vg(:,1), 'r')
hold on
grid on
plot(w, Vg(:,2), 'g')
plot(w, Vg(:,3), 'b')
title('Velocidades de Grupo de los Modos (V_{gj}(\omega))')
xlabel('\omega (rad/s)')
ylabel('V_{gj}(\omega) (m/s)')
legend('V_{g1}(\omega)','V_{g2}(\omega)','V_{g3}(\omega)')
hold off