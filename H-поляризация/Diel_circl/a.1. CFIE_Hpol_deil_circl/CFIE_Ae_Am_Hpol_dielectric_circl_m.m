%%  МИУ

close all;
clear variables;
%% 
% Входные параметры

%  1 - входные данные 
k0 = 2*pi;
% Параметры среды-1
e1 = 2.56;
mu1 = 1;
k1 = k0*sqrt(mu1*e1);
eta1 = 120*pi*sqrt(mu1/e1);

% Параметры среды-2
k2 = k0;
eta2 = 120*pi;

% Параметры программы
a = 5/(2*pi);  % радиус круга
N = 200; % количество точек на круга
phi_i_grad = 0; % угол падения [град]
phi_i_rad = 2*pi/360*phi_i_grad; % угол падения [рад]
gamma = 1.781072417990;

% ========================================================================
%  Построение круга
dphi_grad = 360/N; % угловой шаг на круга [град]
dphi_rad =  2*pi/360*dphi_grad; % угловой шаг на круга [рад]
dx = 2*a*sin(dphi_rad/2);
dx2 = dx/2; 

% круга в полярной системе координат (пусть будет)
% получение векторов углов и вектора ридиусов
phi_circl_rad = zeros(1, N);
phi_circl_grad = zeros(1, N);
a_for_graph = zeros(1, N);
for i = 1 : N
    phi_circl_rad(i) = (i-1)*dphi_rad + dphi_rad/2;
    phi_circl_grad(i) = (i-1)*dphi_grad + dphi_grad/2;
    a_for_graph(i) = a;
end
%% 
% Генерация точек

dphi_grad = 360/N;
x = zeros(1,N);
z = zeros(1,N);

for i = 1 : N
    z(i) = a*sind(dphi_grad*(i-1));
    x(i) = a*cosd(dphi_grad*(i-1));
end

plot(x,z)
%% 
% Посчитаем t для каждого отрезка

for i = 1 : N
    if i == N
        tx(i) = -((x(1)-x(i)))/dx;
        tz(i) = -((z(1)-z(i)))/dx;
    else
        tx(i) = -(x(i+1)-x(i))/dx;
        tz(i) = -(z(i+1)-z(i))/dx;
    end
end
%% 
% Посчитаем n для каждого отрезка

%  Посчитаем n для каждого отрезка
%  нормальный вектор найдем с помощью матрицы повороты
    alpha = -90;
    matrix = [cosd(alpha),-sind(alpha); sind(alpha),cosd(alpha)];

for i = 1 : N
     matrix_temp = [tx(i), tz(i)]*matrix;
     nx(i) = matrix_temp(1);
     nz(i) = matrix_temp(2);
end
%% 
% Координаты центров отрезков

for i = 1 : N 
    if i == N
        x_midle(i) = ((x(1)+x(i)))/2;
        z_midle(i) = ((z(1)+z(i)))/2;
    else
        x_midle(i) = (x(i+1)+x(i))/2;
        z_midle(i) = (z(i+1)+z(i))/2;
    end    
end

plot(x,z,x_midle, z_midle, 'rx');

% визуализируем tx и tz
plot(tx, tz, 'rx');
compass(tx(1),tz(1));

% ========================================================================
% 3 - cам ходя расчета 
% списик для хранения элементов системы 
Z = zeros(2*N,2*N);
Y = zeros(2*N,2*N);
    
for m  = 1 : N
    xm = x_midle(m);
    zm = z_midle(m);
    txm = tx(m);
    tzm = tz(m);
    nxm = nx(m);
    nzm = nz(m);
    
    for n = 1 : N
        xn = x_midle(n);
        zn = z_midle(n);
        txn = tx(n);
        tzn = tz(n);
        nxn = nx(n);
        nzn = nz(n);
            
        % Растояние между началом координат О2 и точкой наблюдений (в ГСК)        
        x_mO2 = xm - xn;
        z_mO2 = zm - zn;
        
        % Положим этот вектор на новую СК
        % Координаты этой точки (в ЛСК)
        umn = x_mO2*txn + z_mO2*tzn;
        wmn = x_mO2*nxn + z_mO2*nzn;
        
        
        % локальные замены 
        xmn = xm - xn;
        zmn = zm - zn;
        r = sqrt(xmn^2 + zmn^2);
        
        umn_plus = umn + dx/2;
        umn_minus = umn - dx/2; 
        r_minus = sqrt(umn_minus^2 + wmn^2);
        r_plus = sqrt(umn_plus^2 + wmn^2);
        rmn = sqrt(umn^2 + wmn^2);
        
        if n == m
              
        % cоставляющие Hy1
        he_y1 =  +1/2;
        hm_y1 = 1i/4*dx*(1 + 1i*2/pi*(log(gamma*k1*dx/4)-1));
        % cоставляющие Hy2  
        he_y2 = -1/2;
        hm_y2 = 1i/4*dx*(1 + 1i*2/pi*(log(gamma*k2*dx/4)-1));
        % cоставляющие E1u
        fe_u12 = 1i/4*dx*(1 + 1i*2/pi*(log(gamma*k2*dx/4)-1));
        fm_u1 = + 1/2; 
        % cоставляющие E2u
        fe_u22 = 1i/4*dx*(1 + 1i*2/pi*(log(gamma*k2*dx/4)-1));
        fm_u2 = - 1/2; 
        % cоставляющие E1w
        fm_w1 = 0;
        % cоставляющие E2w
        fm_w2 = 0;
        
        else
            % cоставляющие Hy1
            he_y1 = - 1i/4*k1*dx*(wmn/rmn*H(1,k1*rmn));
            hm_y1 = 1i/4*dx*H(0,k1*rmn);
            % cоставляющие Hy2  
            he_y2 = - 1i/4*k2*dx*(wmn/rmn*H(1,k2*rmn));
            hm_y2 = 1i/4*dx*H(0,k2*rmn);
            % cоставляющие E1u
            fe_u12 = 1i/4*dx*H(0,k1*rmn);
            fm_u1 =  - 1i/4*k1*dx*(wmn/rmn*H(1,k1*rmn));
            % cоставляющие E2u
            fe_u22 = 1i/4*dx*H(0,k2*rmn);
            fm_u2 = - 1i/4*k2*dx*(wmn/rmn*H(1,k2*rmn)); 
            % cоставляющие E1w
            fm_w1 = + 1i/4*(H(0,k1*r_plus) - H(0,k1*r_minus));
            % cоставляющие E2w
            fm_w2 = + 1i/4*(H(0,k2*r_plus) - H(0,k2*r_minus));
        end
        
      % магнитное поле
        % cоставляющие H1u
        fe_u11 = - 1i/4*k1*(umn_plus/r_plus*H(1,k1*r_plus) - umn_minus/r_minus*H(1,k1*r_minus));
        
        % cоставляющие H2u
        fe_u21 = - 1i/4*k2*(umn_plus/r_plus*H(1,k2*r_plus) - umn_minus/r_minus*H(1,k2*r_minus));
        
        % cоставляющие H1w
        fe_w1 = - 1i/4*k1*(wmn/r_plus*H(1,k1*r_plus) - wmn/r_minus*H(1,k1*r_minus));
        
        % cоставляющие H2w
        fe_w2 = - 1i/4*k2*(wmn/r_plus*H(1,k2*r_plus) - wmn/r_minus*H(1,k2*r_minus));
               
%%
        %% СОБЕРЕМ ВСЕ ВМЕСТЕ
        % cоставляющие Hy1
        H1y = 1i*k1/eta1*hm_y1*eta1 - he_y1;
        
        % cоставляющие Hy2  
        H2y = 1i*k2/eta2*hm_y2*eta2 + he_y2;
        
        % cоставляющие E1u
        E1u = -1i*eta1/k1*(fe_u11 + k1^2*fe_u12) + eta1*fm_u1;
        
        % cоставляющие E2u
        E2u = 1i*eta2/k2*(fe_u21 + k2^2*fe_u22) + eta2*fm_u2;
        
        % cоставляющие E1w
        E1w = -1i*eta1/k1*fe_w1 - eta1*fm_w1;
        
        % cоставляющие E2w
        E2w = 1i*eta2/k2*fe_w2 - eta2*fm_w2;  
        
%%
        % элементы системы       
        H1mn = H1y;
        H2mn = -H2y;
        
        dot_tm_tn = (txm*txn+tzm*tzn);
        dot_tm_nn = (txm*nxn+tzm*nzn);
        
        E1mn = dot_tm_tn*E1u + dot_tm_nn*E1w;
        E2mn = -(dot_tm_tn*E2u + dot_tm_nn*E2w);
        
        % запись в массив для дальнешего решения 
        Z(m,n) = E1mn;
        Z(m,n+N) = E2mn;
        Z(m+N,n) = H1mn;
        Z(m+N,n+N) = H2mn;
    end
end

E = zeros(2*N,1);
% падающее поле 
for m = 1 : N
    xm = x_midle(m);
    zm = z_midle(m);
    txm = tx(m);
    tzm = tz(m);

    % cоставляющие полей
    cx = -1i*k2*cosd(phi_i_grad);
    cz = -1i*k2*sind(phi_i_grad);
    
    Hiy = exp(cx*xm+cz*zm);
    Eix = - eta2*sind(phi_i_grad) * Hiy;
    Eiz = + eta2*cosd(phi_i_grad) * Hiy;
    
    % элементы системы  
    Hi = Hiy;
    Ei = txm*Eix + tzm*Eiz;

    % запись в массив для дальнешего решения 
    E(m) = Ei;
    E(m+N) = Hi;
end   

% Расчитаем токи
% тут зашито сразу 2 тока
I = Z\E;

% отделим один ток от другого
j1 = I(1:N);
j2 = I(N+1:2*N);
%%
% построим эти графики 
plot(phi_circl_grad,  abs(j1), 'r',phi_circl_grad, abs(j2), 'b');

% полярный графки токов
polarplot(phi_circl_rad, abs(j1), 'r',phi_circl_rad, abs(j2), 'b');
%% ЭПР

for p = 1 : 721
    Sum_H = 0;
    phi = (p-1)/2; % угол из формулы для приведения 721 к 360 градусам
    phi_for_graf_DA(p) = (p-1)/2; % угол для графика
    
    cx = - 1i*k2*cosd(phi);
    cz = - 1i*k2*sind(phi);
    

        for n = 1 : N
            txn = tx(n);
            tzn = tz(n);
            xn = x_midle(n);
            zn = z_midle(n);
            
            Sum_H = Sum_H - 1/eta2*eta2*k2/4*dx*j2(n)*exp(cx*xn+cz*zn)...
                    - k2*dx/4*j2(n)*(txn*sind(phi) - tzn*cosd(phi))*exp(cx*xn+cz*zn);    

        end
    
    RCS(p) = 10*log10((4/k2)*Sum_H*conj(Sum_H));
end

% график поля в дальней зоне от угла налюденя (не нормированный)
figure;
plot(phi_for_graf_DA, RCS,'r-');
title('График распределения поля в дальней зоне'); 
xlim([0 360])
xlabel('Угол, град'); 
ylabel('ЭПР');

%% МАГНИТНОЕ Поле на поврехности 
H_surf = zeros(N, 1); 
for m  = 1 : N
    xm = x_midle(m);
    zm = z_midle(m);
    txm = tx(m);
    tzm = tz(m);
    nxm = nx(m);
    nzm = nz(m);
    n_for_graf(m) = m; % для графика
    
    
    H_surf_n = 0;
    for n = 1 : N
        xn = x_midle(n);
        zn = z_midle(n);
        txn = tx(n);
        tzn = tz(n);
        nxn = nx(n);
        nzn = nz(n);
            
        % Растояние между началом координат О2 и точкой наблюдений (в ГСК)        
        x_mO2 = xm - xn;
        z_mO2 = zm - zn;
        
        % Положим этот вектор на новую СК
        % Координаты этой точки (в ЛСК)
        umn = x_mO2*txn + z_mO2*tzn;
        wmn = x_mO2*nxn + z_mO2*nzn;
        
        % локальные замены 
        xmn = xm - xn;
        zmn = zm - zn;
        r = sqrt(xmn^2 + zmn^2);
        
        umn_plus = umn + dx/2;
        umn_minus = umn - dx/2; 
        r_minus = sqrt(umn_minus^2 + wmn^2);
        r_plus = sqrt(umn_plus^2 + wmn^2);
        rmn = sqrt(umn^2 + wmn^2);   
        
        if n == m
            % электрическое поле
            he_y2 = -1/2;
            hm_y2 = 1i/4*dx*(1 + 1i*2/pi*(log(gamma*k2*dx/4)-1));

        else
            % электрчиеское поле
            he_y2 = - 1i/4*k2*dx*(wmn/rmn*H(1,k2*rmn));
            hm_y2 = 1i/4*dx*H(0,k2*rmn);
        end
        H2y = 1i*k2/eta2*hm_y2*eta2 + he_y2;
        H_surf_n = H_surf_n + j2(n) * H2y;
    end
    
    H_surf(m) =  H_surf_n; 
 end
 
figure;
plot(n_for_graf, abs(H_surf), 'r');
title('Электрическое поле не поверхности'); 

% name_file = strcat('I_MIE_circl_diel_Epol_j1_', num2str(N),'.dat');
% f01 = fopen(name_file,'w');
% 
% % за максимум Y примем значения в центре передней стороны падения
% 
% for n=1:N
%     fprintf(f01,' %10.5f %10.5f\n', phi_circl_grad(n), abs(j1(n)));
% end
% 
% fclose(f01);
%%
% name_file = strcat('I_MIE_circl_diel_Epol_j2_', num2str(N),'.dat');
% f02 = fopen(name_file,'w');
% 
% % за максимум Y примем значения в центре передней стороны падения
% 
% for n=1:N
%     fprintf(f02,' %10.5f %10.5f\n', phi_circl_grad(n), abs(j2(n)));
% end
% 
% fclose(f02);
%%
% name_file = strcat('DA_MIE_circl_diel_Epol_j2_', num2str(N),'.dat');
% f03 = fopen(name_file,'w');
% 
% % за максимум Y примем значения в центре передней стороны падения
% 
% for n= 1 : 721
%     fprintf(f03,' %10.5f %10.5f\n', phi_for_graf_DA(n), RCS(n));
% end
% 
% fclose(f03);
%%
function [H_calc] = H(n,x)
    H_calc = besselh(n,1,x); 
end

function [J_calc] = J(n,x)
    J_calc = besselJ(n,x); 
end
    
function [dJ_calc] = dJ(n,x)
    dJ_calc = n/x*besselj(n,x) - besselj(n+1,x);
end

function [dH_calc] = dH(n,x)
    dH_calc = n/x*besselh(n,1,x) - besselh(n+1,1,x);
end