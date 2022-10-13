%+------------------------------------------------------------------------+
%                              Ricardo Vega
%                          SEPI - ESIME Zacatenco
%                  M. en C. en Ingeniería de Telecomunicaciones
%                          Seminario Departamental
%                     Método SPPS para Ecuación de Dispersión
%+------------------------------------------------------------------------+

clear all;
% clc;
close all;

%% Variables a considerar para obtener las soluciones
w = 83E6; %Frecuencia de operación
c = 3E8; %Velocidad de la luz en el espacio libre.
H = 182.09; %Alturas e intervalo de integración
H1 = 255;
H1_H = H1-H;
wp = 2*pi*4.429E6; %Frecuencia de plasma (ionósfera)
ka = w/c; %Número de modo onda en atmósfera 
ki = sqrt(w^2-wp^2)/c; %Número de modo onda en la ionósfera
dim = 5000; %Dimensión vectorial
N = 50; %Número de potencias a considerar
x0 = linspace(0,H,dim/2);
x = linspace(H, H1, dim);
alfa = linspace(ki, ka, dim+2);
alfa = alfa(2:end-1);

%------------------Para el caso 0 <= z <= H -evaluada en H-----------------
u0 = sin(sqrt(ka^2-alfa.^2)'*x0)./sqrt(ka^2-alfa.^2)';   %[alfas x]
u0Der = cos(sqrt(ka^2-alfa.^2)'*x0);

%% Para el caso H <= z <= H1

% Solución de la ecuación homogenea
p = ones(1, dim);
q = zeros(1, dim);
r = -1*ki^2.*ones(1, dim);
lambda = 1;
W = 1;
WDer = 0;

% Construcción de potencias formales
Xvirg = ones(1, dim);
X = ones(1, dim);

% Para "v"
for k = 1: N
    if rem(k, 2) == 1
        Xvirg(k+1, :) = ninteg(Xvirg(k, :)*W^2.*r, H1_H);
        X(k+1, :) = ninteg(X(k, :)./(W^2.*p), H1_H);
    else
        X(k+1, :) = ninteg(X(k, :)*W^2.*r, H1_H);
        Xvirg(k+1, :) = ninteg(Xvirg(k, :)./(W^2.*p), H1_H);
    end
end

% ----------Construcción de soluciones-------------------
v1 = zeros(1, dim);
v2 = zeros(1, dim);
v1Der = zeros(1, dim);
v2Der = zeros(1, dim);

for k = 0: N/2
    v1 = v1 + W*lambda.*Xvirg(2*k+1, :);
    if k > 0
        v1Der = v1Der + lambda*Xvirg(2*k, :)./(W*p);
    end
end                                       

for k = 0: N/2-1
    v2 = v2 + lambda*W.*X(2*k+2, :);
    v2Der = v2Der + lambda*X(2*k+1, :)./(W*p);
end

% Para "u"
v = v1 + 1i*v2;
vDer = v1Der + 1i*v2Der;

p = ones(1, dim);
q = ki^2.*ones(1, dim);
r = ones(1, dim);

for k = 1: N 
    if rem(k, 2) == 1
        Xvirg(k+1, :) = ninteg(Xvirg(k, :).*v.^2.*r, H1_H);
        X(k+1, :) = ninteg(X(k, :)./(v.^2.*p), H1_H);
    else
        X(k+1, :) = ninteg(X(k, :).*v.^2.*r, H1_H);
        Xvirg(k+1, :) = ninteg(Xvirg(k, :)./(v.^2.*p), H);
    end
end

%-------------------------------------------------------------------------------------------------------
%Construcción de la solución en H <= x <= H1

u1 = zeros(N/2+1,dim);  %[potencias x]
u2 = zeros(N/2+1,dim);
u1Der = zeros(N/2+1,dim);
u2Der = zeros(N/2+1,dim);

for cont = 0:N/2
    u1(cont+1,:) = v.*Xvirg(1+2*cont,:);
    if cont > 0
        u1Der(cont+1,:) = Xvirg(2*cont,:)./(v.*p);
    end
end

u1Der = u1.*vDer./v + u1Der;
u1 = flip(u1);
u1Der = flip(u1Der);

for cont = 0:N/2-1
    u2(cont+1,:) = v.*X(2+2*cont,:);
    u2Der(cont+1,:) = X(1+2*cont,:)./(v.*p);
end

u2Der = u2.*vDer./v + u2Der;
u2 = flip(u2);
u2Der = flip(u2Der);

%---------------Prueba para x, con alfa constante--------------------
NoAlfa = 4999;     %se eligen valores entre 1 y 5000

for cont = 1:dim
    u1x(cont) = polyval(u1(:,cont),alfa(NoAlfa).^2); %Funciones de x
    u2x(cont) = polyval(u2(:,cont),alfa(NoAlfa).^2);
    u1Derx(cont) = polyval(u1Der(:,cont),alfa(NoAlfa).^2);
    u2Derx(cont) = polyval(u2Der(:,cont),alfa(NoAlfa).^2);
end

C1 = u0(NoAlfa,end)/u1x(1);
C2 = u0Der(NoAlfa,end)/u2Derx(1);
uiAC = C1*real(u1x) + C2*real(u2x);
uiDerAC = C1*real(u1Derx) + C2*real(u2Derx);

% figure
% plot(x0,u0(NoAlfa,:),'b')
% hold on
% plot(x0,u0Der(NoAlfa,:),'r')
% plot(x,uiAC,'m');
% plot(x,uiDerAC,'g');
% grid on
% title('Solutions for a particular \alpha')
% legend('u_{0}', 'u_{0}\prime', 'u_{i}','u_{i}\prime')
% xlabel('x [km]')


%----------------Construcción de fii y fiiDer---------------------------
u1Val = zeros(dim,dim);    %[alfas x]
u2Val = zeros(dim,dim);
u1DerVal = zeros(dim,dim);
u2DerVal = zeros(dim,dim);

for cont = 1:dim        %Evaluación para cada punto en x (columnas)
    u1Val(:,cont) = polyval(u1(:,cont),alfa.^2); %renglones con valores para alfa
    u2Val(:,cont) = polyval(u2(:,cont),alfa.^2);
    u1DerVal(:,cont) = polyval(u1Der(:,cont),alfa.^2);
    u2DerVal(:,cont) = polyval(u2Der(:,cont),alfa.^2);
end

C1 = u0(:,end)./u1Val(:,1);
C2 = u0Der(:,end)./u2DerVal(:,1);

ui = real(u1Val).*C1 + real(u2Val).*C2;
uiDer = real(u1DerVal).*C1 + real(u2DerVal).*C2;

%--Figura para verificar unión de soluciones(alfa) para altura H--

% figure
% plot(alfa.^2,fi0(:,end),'b')
% hold on
% plot(alfa.^2,fi0Der(:,end),'r')
% plot(alfa.^2,fii(:,1),'c--')
% plot(alfa.^2,fiiDer(:,1),'y--')
% grid on
% title('Evaluation at z = H')
% xlabel('\alpha^{2}')
% legend('\phi_{0}','\phi_{0}\prime','\phi_{i}','\phi_{i}\prime')

%-------------------Ecuación de Dispersión-----------------------------
EcDisp = ui(:,1) + uiDer(:,1)./sqrt(alfa.^2-ki^2)';
EcDispSpline = spline(alfa.^2,EcDisp);
Eigenvalores = fnzeros(EcDispSpline); %localiza los ceros del spline
NMO = sqrt(Eigenvalores(1,:))

figure
plot(alfa, EcDisp)
hold on
% plot(alfa, ppval(EcDispSpline,alfa.^2),'c--')
stem(NMO, ones(1, length(NMO)));
grid on
title('\Lambda(\omega,\alpha)')
xlabel('\alpha')
legend('\Lambda', '\alpha_{n}')