classdef D_class_MIE_EFIE_sqr_H_pol < P_class_MIE_sqr
%% """Класс для MIE EFIE E-pol sqr"""
% 02-03-21

% Пример использования:
    % as = 1;
    % phi_grad = 60;
    % N= 100;
    %     
    % calc_1 = D_class_MIE_EFIE_sqr_H_pol(as, phi_grad);
    % calc_1.calculate_main_sqr(N)
    % calc_1.calculate_I(N).write_I_to_file();
    % calc_1.calculate_DA().write_DA_to_file();
    % I = calc_1.get_I();
    % DA = calc_1.get_DA();
    
%% *** СВОЙСТВА КЛАССА ***
%%  Уже расчитанные значения    
    %%  Уже расчитанные значения
    properties 
        % пусто пока что - все есть в родительском классе
    end   
    
%% *** МЕТОДЫ***
    methods  
%% КОНСТРУКТОР - НАСЛЕДОВАНИЕ
        % конструктор
        function obj = D_class_MIE_EFIE_sqr_H_pol(as, phi_grad)
            obj@P_class_MIE_sqr(as, phi_grad);
            obj.type_MIE = "EFIE";
            obj.type_pol = "Hpol";
            obj.name_of_class = strcat("class\_EFIE\_Hpol\_sqr");
        end
%% Генерация D0        

        % Генерация D0      
        function D0 = ganerate_D0(obj)
            % Генерация D0
            % ВЫТАЩИМ
            dx = obj.dx;
            k = obj.k;
            gamma = obj.gamma;
            
            % ввели значение
            kd4 = k*dx/4;
            D0 = k*dx*(1+(1i*2/pi)*(log(gamma*kd4)-1));
        end     
        
%% РАСЧЕТ ТОКОВ       
        % расчет токов для данного числа источников и точек на элипсе
        function obj = calculate_I(obj, N)
            
            % Изучим один цикл
            % генерация точек для эллипса
            obj.calculate_main_sqr(N);
            % достать параметры от сюда
            % достать параметры от сюда
            X = obj.sqr_main.X;
            Z = obj.sqr_main.Z;
            tx = obj.sqr_main.tx;
            tz = obj.sqr_main.tz;
            nx = obj.sqr_main.nx; 
            nz = obj.sqr_main.nz; 
            dx = obj.dx;
            k = obj.k;
            phi_grad = obj.phi_grad;
            gamma = obj.gamma;
            
            % построим эллипс
            obj.print_sqr()
            
            % cгенерируем DO для расчета 
            D0 = ganerate_D0(obj);
            
            % определим замены 
            kdx = k*dx;
            delta = dx/2;
            
            % Изучим один цикл
            count = 0;
            for m = 1 : 4*N
            %   Координаты в исходной системе координат (в ГСК) 
            % - т.к. точки находятся в центре отрезков разбиения
                xm = X(m);
                zm = Z(m);
                txm = tx(m);
                tzm = tz(m);
                nxm = nx(m);
                nzm = nz(m);

                for n = 1 : 4*N
                %   Координаты в исходной системе координат (в ГСК)
                    xn = X(n);
                    zn = Z(n);
                    txn = tx(n);
                    tzn = tz(n);
                    nxn = nx(n);
                    nzn = nz(n);

                %   Центр новой системы кооридны - начало отрезка (в ГСК)
                %   Давай возьмем середину отрезка
                    x_begin_O2 = xn;
                    z_begin_O2 = zn;

                %   Растояние между началом координат О2 и точкой наблюдений (в ГСК)
                    x_mO2 = xm - x_begin_O2;
                    z_mO2 =  zm - z_begin_O2;

            %       Положим этот вектор на новую СК
            %       Координаты этой точки (в ЛСК)
                    xmO2 = x_mO2*txn + z_mO2*tzn;
                    zmO2 = x_mO2*nxn + z_mO2*nzn;

            %       Т.к. взяли на начало середину отрезка то тут будет
                    % координаты в локальных координатах
%                     umn = xmO2;
%                     wmn = zmO2;                    
                    xmnO2 = xmO2;
                    zmnO2 = zmO2;

            %       Растояние между m и n (в ГСК)
                    xmn = xm - xn;
                    zmn = zm - zn;
                    r = sqrt(xmn^2 + zmn^2);

            %       Расчет матричных элементов (в ЛСК)
                    if n == m %r < 0.001*dx 
                        Zm = D0;

                    elseif abs(n - m) < 1.5 %r <= 2*dx
                        xmn_plus_delta = (xmnO2) + delta;
                        xmn_minus_delta = (xmnO2) - delta;
                        r_plus_delta = sqrt(xmn_plus_delta^2 + zmnO2^2);
                        r_minus_delta = sqrt(xmn_minus_delta^2 + zmnO2^2);
                        zmn = zmnO2;
                        Zm = kdx+k*2*1i/pi*(dx*log(gamma*k/2)-dx*zmn*(atan(xmn_plus_delta/zmn)- atan(xmn_minus_delta/zmn)) + 0.5*xmn_plus_delta*log(r_plus_delta^2)  -  0.5* xmn_minus_delta*log(r_minus_delta^2)- dx);
                    else
                        f2 = besselh(0,1,k*r);
                        Zm = k*dx*f2;
                    end
                    
                        r_x_plus = sqrt((xmnO2+delta)^2+zmnO2^2);
                        r_x_minus = sqrt((xmnO2-delta)^2+zmnO2^2);

                        Esx= (xmnO2-delta)/r_x_minus*besselh(1, 1, k*r_x_minus) - (xmnO2+delta)/r_x_plus*besselh(1, 1, k*r_x_plus) + Zm;
                        Esz = (zmnO2)/r_x_minus*besselh(1, 1, k*r_x_minus) - (zmnO2)/r_x_plus*besselh(1, 1, k*r_x_plus);        

                        tx_betven_vectors = txm*txn + tzm*tzn;
                        tz_betven_vectors = txm*nxn + tzm*nzn;

                        Es(m,n) = Esx*tx_betven_vectors + Esz*tz_betven_vectors;
                end
            end

            % Расчет падающего поля
            for m=1: 4*N
                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);
                % Определим координаты
                xm = X(m);
                zm = Z(m);
                txm = tx(m);
                tzm = tz(m);

            %   Итоговое падающее поел - это суперпозиция полей
                Eix = +4*1i*sind(phi_grad)*exp(cx*xm+cz*zm);% падающее поле X
                Eiz = -4*1i*cosd(phi_grad)*exp(cx*xm+cz*zm);% падающее поле Z
                Ei(m,1) = Eix*txm + Eiz*tzm; 
            end

            % Посчитаем токи и построим график
            % расчеток токов ВИ 
            I=Es\Ei;
            
            % номера ВИ для графика 
            n_for_graf = zeros(1, N);
            for i = 1 : 4*N
                 n_for_graf(i) = i;
            end
            
            % ось X в относительных единицах
            figure;
            plot(n_for_graf*dx, abs(I), 'r');
            title('График распределения тока (ток около 2)  - (4)'); 
            xlabel('Периметр'); 
            ylabel('Плотность тока');
            legend(obj.name_of_class);
            
            % ЗАПИШЕМ необходимые параметры из СВ-В КЛАССА
            obj.I_real = I; 
            obj.n_for_graf = n_for_graf;
            obj.Zmn = Es;
            obj.Ei = Ei;
           
        end
        
%% РАСЧЕТ ЭПР       
        % расчитать ЭПР
        function obj = calculate_DA(obj)
            
            % ВЫТАЩИМ необходимые параметры из СВ-В КЛАССА 
            k = obj.k;
            I_real = obj.I_real;
            N = obj.N;
            tx = obj.sqr_main.tx;
            tz = obj.sqr_main.tz;
            X = obj.sqr_main.X;
            Z = obj.sqr_main.Z;
            dx = obj.dx;
            
            for p = 1 : 721
                Sum_H = 0;
                phi = (p-1)/2; % угол из формулы для приведения 721 к 360 градусам
                phi_for_graf_DA(p) = (p-1)/2; % угол для графика

                cx = - 1i*k*cosd(phi);
                cz = - 1i*k*sind(phi);


                    for n = 1 : 4*N
                        txn = tx(n);
                        tzn = tz(n);
                        xn = X(n);
                        zn = Z(n);
            
                        % ДЕЛИТЬ ЛИ НА 120 ПИ *(120*pi)
                        Sum_H = Sum_H + 1/4*I_real(n)*k*dx*(txn*sind(phi) - tzn*cosd(phi))*exp(cx*xn+cz*zn);

                    end
                    
            Sum_H_mass(p) = Sum_H;
            RCS(p) = 10*log10((2/pi)*Sum_H*conj(Sum_H));

            end

            % график поля в дальней зоне от угла налюденя (не нормированный)
            figure;
            plot(phi_for_graf_DA ,RCS,'r-');
            title('График распределения поля в дальней зоне'); 
            xlim([0 360])
            xlabel('Угол, град'); 
            ylabel('ЭПР');
            legend(obj.name_of_class);
            
            % ЗАПИШЕМ
            obj.DA = RCS; 
            obj.phi_for_graf_DA = phi_for_graf_DA;
            obj.Sum_H_mass = Sum_H_mass;
            
        end
    end    
end