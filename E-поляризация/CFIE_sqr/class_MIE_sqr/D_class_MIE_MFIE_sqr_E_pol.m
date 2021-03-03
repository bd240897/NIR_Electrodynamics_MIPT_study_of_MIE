classdef D_class_MIE_MFIE_sqr_E_pol < P_class_MIE_sqr
%% """Класс для MIE EFIE E-pol sqr"""
% 02-03-21

% Пример использования:
    % as = 1;
    % phi_grad = 60;
    % N= 100;
    %     
    % calc_1 = D_class_MIE_MFIE_sqr_E_pol(as, phi_grad);
    % calc_1.calculate_main_sqr(N)
    % calc_1.calculate_I(N).write_I_to_file();
    % calc_1.calculate_DA().write_DA_to_file();
    % I = calc_1.get_I();
    % DA = calc_1.get_DA();
    
%% *** СВОЙСТВА КЛАССА ***
%%  Уже расчитанные значения    
    %%  Уже расчитанные значения
    properties 

    end   
    
%% *** МЕТОДЫ***
    methods  
%% КОНСТРУКТОР - НАСЛЕДОВАНИЕ
        % конструктор
        function obj = D_class_MIE_MFIE_sqr_E_pol(as, phi_grad)
            obj@P_class_MIE_sqr(as, phi_grad);
            obj.type_MIE = string("MFIE");
            obj.type_pol = string("Epol");
            obj.name_of_class = strcat("class\_MFIE\_Epol\_sqr");
        end
%% РАСЧЕТ ТОКОВ       
        % расчет токов для данного числа источников и точек на элипсе
        function obj = calculate_I(obj, N)

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
            type_MIE = obj.type_MIE;
            type_pol = obj.type_pol;
            
            % построим эллипс
            obj.print_sqr()
                               
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
                    umn = xmO2;
                    wmn = zmO2;                    

            %       Растояние между m и n (в ГСК)
                    xmn = xm - xn;
                    zmn = zm - zn;
                    r = sqrt(xmn^2 + zmn^2);

            %       Расчет матричных элементов (в ЛСК)
                    Hus = 0;
                    Hws = 0;
                    
                    if n == m %r < 0.001*dx 

                        % u - cоставляющая 
                        Hus = 0.5;

                        % w - cоставляющая
                        Hws = 0;

                    elseif abs(n - m) < 1.5 %r <= 2*dx
                        % u - cоставляющая                     
                        umn_plus_dx2 = umn + dx/2;
                        umn_minus_dx2 = umn - dx/2;
                        Hus = (1i/8*wmn*k^2*dx + 1/(2*pi)*(atan(umn_plus_dx2/wmn) - atan(umn_minus_dx2/wmn)));

                        % w - cоставляющая
                        r_plus_dx2 = sqrt((umn + dx/2)^2 + wmn^2);
                        r_minus_dx2 = sqrt((umn - dx/2)^2 + wmn^2);
                        Hws = 1i/4*(besselh(0,1,k*r_plus_dx2) - besselh(0,1,k*r_minus_dx2));

                    else
                        % u - cоставляющая
                        rmn = sqrt(umn^2 + wmn^2);    
                        Hus = 1i/4*k*dx*(wmn/rmn*besselh(1,1,k*rmn));  

                        % w - cоставляющая
                        r_plus_dx2 = sqrt((umn + dx/2)^2 + wmn^2);
                        r_minus_dx2 = sqrt((umn - dx/2)^2 + wmn^2);
                        Hws = 1i/4*(besselh(0,1,k*r_plus_dx2) - besselh(0,1,k*r_minus_dx2));
                    end

            %       символ кронкера
                    if m == n
                        kronker = 1;
                    else 
                        kronker = 0;
                    end

                    % расчет матричнх элементов 
                    Zmn(m,n) = kronker + Hus*(nxm*tzn - nzm*txn) + Hws*(nxm*nzn - nzm*nxn);

                end
            end

            % Расчет падающего поля
            for m= 1 : 4*N
                nxm = nx(m);
                nzm = nz(m);

                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);

                % Определим координаты
                xm = X(m);
                zm = Z(m);

            %   Итоговое падающее поел - это суперпозиция полей
            %     etta = 120*pi;
                Hxi = 1*sind(phi_grad)*exp(cx*xm+cz*zm);
                Hzi = -1*cosd(phi_grad)*exp(cx*xm+cz*zm);

                Hi(m,1) = -nxm*Hzi + nzm*Hxi;% падающее поле X
            end

            % Посчитаем токи и построим график
            % расчеток токов ВИ 
            I=Zmn\Hi;
            
            % номера ВИ для графика 
            n_for_graf = zeros(1, N);
            for i = 1 : 4*N
                 n_for_graf(i) = i;
            end
            
            % ось X в относительных единицах
            figure;
            plot(n_for_graf*dx, abs(I), 'r');
            

            
            
            title("График распределения тока (ток около 2)  - (4)"); 
            xlabel("Периметр"); 
            ylabel("Плотность тока");
            legend(obj.name_of_class);
            
            % ЗАПИШЕМ необходимые параметры из СВ-В КЛАССА
            obj.I_real = I; 
            obj.n_for_graf = n_for_graf;
            obj.Zmn = Zmn;
            obj.Hi = Hi;
            
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
            type_MIE = obj.type_MIE;
            type_pol = obj.type_pol;
            
            for p = 1 : 721
                Sum_E = 0;
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
                        Sum_E = Sum_E + k/4*dx*I_real(n)*exp(cx*xn+cz*zn);      

                    end
                    
            Sum_E_mass(p) = Sum_E;
            RCS(p) = 10*log10((2/pi)*Sum_E*conj(Sum_E));

            end

            % график поля в дальней зоне от угла налюденя (не нормированный)
            figure;
            plot(phi_for_graf_DA ,RCS,'r-');
            title("График распределения поля в дальней зоне"); 
            xlim([0 360])
            xlabel("Угол, град"); 
            ylabel("ЭПР");
            legend(obj.name_of_class);
            
            % ЗАПИШЕМ
            obj.DA = RCS; 
            obj.phi_for_graf_DA = phi_for_graf_DA;
            obj.Sum_E_mass = Sum_E_mass;
           
        end
    end    
end