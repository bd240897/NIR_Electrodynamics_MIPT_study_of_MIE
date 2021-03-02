classdef class_MIE_MFIE_sqr_E_pol < handle
%% """Класс для MIE EFIE E-pol sqr"""
% 02-03-21

% Пример использования:
    % as = 1;
    % phi_grad = 60;
    % N= 100;
    %     
    % calc_1 = class_MIE_MFIE_sqr_E_pol(as, phi_grad);
    % calc_1.calculate_main_sqr(N)
    % calc_1.calculate_I(N).write_I_to_file();
    % calc_1.calculate_DA().write_DA_to_file();
    % I = calc_1.get_I();
    % DA = calc_1.get_DA();
    
    %%  Уже расчитанные значения
    properties 
        sqr_main = struct();  
        I_real = [];
        DA = [];
        % то что будем дсотавать из вне 
        Sum_E_mass = [];
        Zmn = [];
        Hi = [];
    end
%% Параметры элипса и расчета
    properties
        as = NaN;
        phi_grad = NaN;
        N = NaN;
        delta = [];
        dx = [];
        n_for_graf = NaN;
        phi_for_graf_DA = NaN;
    end
    %% Параметры вспомогательного элипса
    properties
        
    end
%%  Константы   
    properties
        k = 2*pi;
        gamma = 1.781072417990;     
    end
    
 %%   
    methods  
%%
        % конструктор
        function obj = class_MIE_MFIE_sqr_E_pol(as, phi_grad)
            obj.as = as;
            obj.phi_grad = phi_grad;
        end
%%        
        % расчет основного эллипса
        % выбор кол-ва точек на элипсе 
        function obj  = calculate_main_sqr(obj, N)
            % ЗАПИШЕМ в свойства N и бдудем ег оиспользовать везде далее
            obj.N = N;
            
            % вытащи необходимые параметры из СВ-В КЛАССА
            % полуоси элипса
            as = obj.as; %cтороная квадрата
            k  = obj.k;
            
            % расчет dx
            dx = as/N;
            
            % Генерация точек
            x = zeros(1,N);
            z = zeros(1,N);

            % Генерация точек
            % Расчет кооодинат центров отрезков разбиения (без функции)
            % генератора координат для первого участка
            for n=1:N
                x1(n) = (n-0.5)*dx;
                z1(n) = 0;
                x2(n) = as;
                z2(n) = (n-0.5)*dx;
                x3(n) = (n-0.5)*dx;
                z3(n) = as;
                x4(n) = 0;
                z4(n) = (n-0.5)*dx;  
            end

            % координта на поверхности цилиндра
            X = [x1 x2 fliplr(x3) x4];
            Z = [z1 z2 z3 fliplr(z4)];

%             plot(X,Z)
%             xlim([-1,3])
%             ylim([-1,3])

            % Расчет касательных к отрезкам разибения
            % касательные * до копирования
            tx_one = [1 0 -1 0];
            tz_one = [0 1 0 -1];

            %  общие касательный
            tz = [repmat(tz_one(1),1,N) repmat(tz_one(2),1,N) repmat(tz_one(3),1,N) repmat(tz_one(4),1,N)];
            tx = [repmat(tx_one(1),1,N) repmat(tx_one(2),1,N) repmat(tx_one(3),1,N) repmat(tx_one(4),1,N)];
            
            % Посчитаем n для каждого отрезка
            %   нормальный вектор найдем с помощью матрицы повороты
                a = -90;
                matrix = [cosd(a),-sind(a); sind(a),cosd(a)];

            for i = 1 : 4*N
                 matrix_temp = [tx(i), tz(i)]*matrix;
                 nx(i) = matrix_temp(1);
                 nz(i) = matrix_temp(2);
            end
            
            % введем что такое дельта
            delta = dx/2;
                       
            % ЗАПИШЕМ необходимые параметры из СВ-В КЛАССА          
            obj.sqr_main.X = X;
            obj.sqr_main.Z = Z;
            obj.sqr_main.tx = tx;
            obj.sqr_main.tz = tz;
            obj.sqr_main.nx = nx;
            obj.sqr_main.nz = nz;
            obj.dx = dx;
            obj.delta = delta;
        end
%% 
        % функция для построения графика эллипса 
        function print_sqr(obj) 
            % ВЫТАЩИМ
            X = obj.sqr_main.X;
            Z = obj.sqr_main.Z;
            X = obj.sqr_main.X;
            Z = obj.sqr_main.Z;
            
            figure;
            plot(X, Z, 'b.');
            title('Основной квадрат и его проксимация (MAS EFIE Epol)'); 

        end
%% 
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
%%     
        % расчет токов для данного числа источников и точек на элипсе
        function obj = calculate_I(obj, N)

            % генерация точек для эллипса
            sqr_main  = obj.calculate_main_sqr(N);
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
            for m=1 : 4*N
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
            title('График распределения тока (ток около 2)  - (4)'); 
            xlabel('Периметр'); 
            ylabel('Плотность тока');
            
            % ЗАПИШЕМ необходимые параметры из СВ-В КЛАССА
            obj.I_real = I; 
            obj.n_for_graf = n_for_graf;
            obj.Zmn = Zmn;
            obj.Hi = Hi;
            
        end
%%        
        % записать токи в файл
        function obj = write_I_to_file(obj)
            
            % ВЫТАЩИМ 
            I_real = obj.I_real;
            N = obj.N;
            n_for_graf = obj.n_for_graf;
            
            name_file = strcat('I_HFIE_MIE_Epol_sqrt_', num2str(N),'.dat');
            
            f01 = fopen(name_file,'w');
            for n=1:N
                fprintf(f01,' %10.5f %10.5f\n', n_for_graf(n)*dx, abs(I_real(n)));
            end         
            fclose(f01);         
        end
%%        
        % расчитать ЭПР
        function obj= calculate_DA(obj)
            
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
            title('График распределения поля в дальней зоне (MIE MFIE Epol)'); 
            xlim([0 360])
            xlabel('Угол, град'); 
            ylabel('ЭПР');
            
            % ЗАПИШЕМ
            obj.DA = RCS; 
            obj.phi_for_graf_DA = phi_for_graf_DA;
            obj.Sum_E_mass = Sum_E_mass;
        end
 %%       
        % взаписать токи в файл
        function obj= calculate_DE(obj)
        end
 %%       
        % расчитать невязку
        function obj =  write_DA_to_file(obj)
            
            % ВЫТАЩИМ
            phi_for_graf_DA = obj.phi_for_graf_DA;
            DA = obj.DA;
            
            f02 = fopen('DA_MFIE_MIE_Epol_sqrt.dat','w');
            for p = 1 : 721
                fprintf(f02, ' %10.5f %10.5f\n', phi_for_graf_DA(p), DA(p));
            end
            fclose(f02);
        end
                %%
        % получить ток -  GETer
        function I = get_I(obj)
            I = obj.I_real;
        end
        % получить ЭПР -  GETer
        function DA = get_DA(obj)
            DA = obj.DA;
        end
    end
end