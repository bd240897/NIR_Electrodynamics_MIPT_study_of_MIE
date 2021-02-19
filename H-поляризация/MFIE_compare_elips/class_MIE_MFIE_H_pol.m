classdef class_MIE_MFIE_H_pol < handle
%% """Класс для MIE MFIE H-pol"""
% 19-02-21

% Пример использования:
    % a = 2;
    % b = 1;
    % phi_grad = 60;
    % N= 100;
    %     
    % calc_1 = class_MIE_MFIE_H_pol(a,b,phi_grad);
    % calc_1.calculate_I(N).write_I_to_file();
    % calc_1.calculate_DA().write_DA_to_file();
    % I = calc_1.get_I();
    % DA = calc_1.get_DA();
    
    %%  Уже расчитанные значения
    properties 
        ellips_main = struct();  
        I = [];
        Sum_H_mass = [];
        DA = [];
    end
%% Параметры элипса и расчета
    properties
        a = NaN;
        b = NaN;
        phi_grad = NaN;
        N = NaN;
        delta = [];
        dx = [];
        phi_for_graf_real = NaN;
        phi_for_graf_DA = NaN;
    end
%%  Константы   
    properties
        k = 2*pi;
        gamma = 1.781072417990;;        
    end
    
 %%   
    methods  
%%
        % конструктор
        function obj = class_MIE_MFIE_H_pol(a,b,phi_grad)
            obj.a = a;
            obj.b = b;
            obj.phi_grad = phi_grad;
        end
%%        
        % расчет основного эллипса
        % выбор кол-ва точек на элипсе 
        function obj  = calculate_main_ellips(obj, N)
            % ЗАПИШЕМ в свойства N и бдудем ег оиспользовать везде далее
            obj.N = N;
            
            % вытащи необходимые параметры из СВ-В КЛАССА
            % полуоси элипса
            a = obj.a; %большшая полуось - ОX
            b = obj.b; %меньшая полуось - ОZ
            k  = obj.k;

            % Генерация точек
            delta_phi = 360/N;
            x = zeros(1,N);
            z = zeros(1,N);

            for i = 1 : N
                z(i) = b*sind(delta_phi*(i-1));
                x(i) = a*cosd(delta_phi*(i-1));
            end

%             plot(x,z)
%             title("Основной эллипс (MIE MFIE Hpol)");
            
            % Расчет dx для каждого участка
            for i = 1 : N
                if i == N
                    d_x = (x(1)-x(i));
                    d_z = (z(1)-z(i));
                else
                    d_x = (x(i+1)-x(i));
                    d_z = (z(i+1)-z(i));
                end
                dx(i) = sqrt(d_x^2+d_z^2);    
            end


            % Посчитаем t для каждого отрезка
            for i = 1 : N
                if i == N
                    tx(i) = -((x(1)-x(i)))/dx(i);
                    tz(i) = -((z(1)-z(i)))/dx(i);
                else
                    tx(i) = -(x(i+1)-x(i))/dx(i);
                    tz(i) = -(z(i+1)-z(i))/dx(i);
                end
            end
            
            % Посчитаем n для каждого отрезка
            %   нормальный вектор найдем с помощью матрицы повороты
                alpha = -90;
                matrix = [cosd(alpha),-sind(alpha); sind(alpha),cosd(alpha)];

            for i = 1 : N
                 matrix_temp = [tx(i), tz(i)]*matrix;
                 nx(i) = matrix_temp(1);
                 nz(i) = matrix_temp(2);
            end

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
            
%             figure;
%             plot(x,z,x_midle, z_midle, 'rx')
%             title('Эллипс и его проксимация (MIE MFIE Hpol)'); 
            
            % визуализируем tx и tz
%             figure;
%             plot(tx, tz, 'rx')
% %             compass(tx(1),tz(1))
%             title('Касательные для элипса (MIE MFIE Hpol)'); 
            
            delta = dx/2;
                       
            % ЗАПИШЕМ необходимые параметры из СВ-В КЛАССА          
            obj.ellips_main.x = x;
            obj.ellips_main.z = z;
            obj.ellips_main.x_midle = x_midle;
            obj.ellips_main.z_midle = z_midle;
            obj.ellips_main.tx = tx;
            obj.ellips_main.tz = tz;
            obj.ellips_main.nx = nx;
            obj.ellips_main.nz = nz;
            obj.dx = dx;
            obj.delta = delta;
        end
%% 
        % функция для построения графика эллипса 
        function print_ellips(obj) 
            % ВЫТАЩИМ
            x = obj.ellips_main.x;
            z = obj.ellips_main.z;
            x_midle = obj.ellips_main.x_midle;
            z_midle = obj.ellips_main.z_midle;
            
            figure;
            plot(x, z, 'r', x_midle, z_midle, 'b.');
            title('Основной эллипс и его проксимация (MAS MFIE Hpol)'); 

        end
%%     
        % расчет токов для данного числа источников и точек на элипсе
        function obj = calculate_I(obj, N)
            % генерация точек для эллипса
            ellips_main  = obj.calculate_main_ellips(N);

            % достать параметры от сюда
            x_midle = obj.ellips_main.x_midle;
            z_midle = obj.ellips_main.z_midle;
            tx = obj.ellips_main.tx;
            tz = obj.ellips_main.tz;
            nx = obj.ellips_main.nx; 
            nz = obj.ellips_main.nz; 
            dx = obj.dx;
            k = obj.k;
            phi_grad = obj.phi_grad;
            
            % построим эллипс
            obj.print_ellips()
            
            % Изучим один цикл
            count = 0;
            for m = 1 : N
            %   Координаты в исходной системе координат (в ГСК) 
            % - т.к. точки находятся в центре отрезков разбиения
                xm = x_midle(m);
                zm = z_midle(m);
                txm = tx(m);
                tzm = tz(m);
                nxm = nx(m);
                nzm = nz(m);

                for n = 1 : N
                %   Координаты в исходной системе координат (в ГСК)
                    xn = x_midle(n);
                    zn = z_midle(n);
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

                    cur_dx = dx(n); % current dx

                    if n == m %r < 0.001*dx 
                        % y - cоставляющая 
                        Hys = -0.5;

                    elseif abs(n - m) < 1.5 %r <= 2*dx
                        % y - cоставляющая                     
                        umn_plus_dx2 = umn + cur_dx/2;
                        umn_minus_dx2 = umn - cur_dx/2;
                        Hys = -(1i/8*wmn*k^2*cur_dx + 1/(2*pi)*(atan(umn_plus_dx2/wmn) - atan(umn_minus_dx2/wmn)));

                    else 
                        % y - cоставляющая
                        rmn = sqrt(umn^2 + wmn^2);    
                        Hys = -1i/4*k*cur_dx*(wmn/rmn*besselh(1,1,k*rmn));  

                    end

            %       символ кронкера
                    if m == n
                        kronker = 1;
                    else 
                        kronker = 0;
                    end

                    % расчет матричнх элементов 
                    Zmn(m,n) = kronker + Hys;
                end
            end

            % Расчет падающего поля
            for m=1:N
                nxm = nx(m);
                nzm = nz(m);

                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);

                % Определим координаты
                xm = x_midle(m);
                zm = z_midle(m);

            %   Итоговое падающее поле
                Hi(m,1) = exp(cx*xm+cz*zm); 
            end

            % Посчитаем токи и построим график
            % расчеток токов ВИ 
            I=Zmn\Hi;
            
            for i = 1 : N
                z_m = z_midle(i);
                x_m = x_midle(i);

                if (x_m <= 0 && z_m < 0) || (x_m >= 0 && z_m < 0)
                    phi_for_graf_real(i) = atan2(z_m, x_m)*180/pi+360;
                else
                    phi_for_graf_real(i) = atan2(z_m, x_m)*180/pi;
                end
            end

            figure;
            fig = plot(phi_for_graf_real, abs(I), 'r');
            title('График распределения тока в ВИ (MIE MFIE Hpol)'); 
            xlabel('Угол'); 
            xlim([0 360])
            ylabel('Плотность тока');
            
            % ЗАПИШЕМ необходимые параметры из СВ-В КЛАССА
            obj.I = I; 
            obj.phi_for_graf_real = phi_for_graf_real;
            
        end
%%        
        % записать токи в файл
        function obj = write_I_to_file(obj)
            
            % ВЫТАЩИМ 
            I = obj.I;
            N = obj.N;
            phi_for_graf_real = obj.phi_for_graf_real;
            
            name_file = strcat('I_MFIE_MIE_Hpol_elips_', num2str(N),'.dat');
            
            f01 = fopen(name_file,'w');
            for n=1:N
                fprintf(f01,' %10.5f %10.5f\n', phi_for_graf_real(n), abs(I(n))*1000);
            end         
            fclose(f01);         
        end
%%        
        % расчитать ЭПР
        function obj= calculate_DA(obj)
            
            % ВЫТАЩИМ необходимые параметры из СВ-В КЛАССА 
            k = obj.k;
            I = obj.I;
            xAS = obj.ellips_main.x_midle;
            zAS = obj.ellips_main.z_midle;
            N = obj.N;
            tx = obj.ellips_main.tx;
            tz = obj.ellips_main.tz;
            x_midle = obj.ellips_main.x_midle;
            z_midle = obj.ellips_main.z_midle;
            dx = obj.dx;
            
            for p = 1 : 721
                Sum_H = 0;
                phi = (p-1)/2; % угол из формулы для приведения 721 к 360 градусам
                phi_for_graf_DA(p) = (p-1)/2; % угол для графика

                cx = - 1i*k*cosd(phi);
                cz = - 1i*k*sind(phi);

                for n = 1 : N
                   txn = tx(n);
                   tzn = tz(n);
                   xn = x_midle(n);
                   zn = z_midle(n);

                   Sum_H = Sum_H + 1/4*I(n)*k*dx(n)*(txn*sind(phi) - tzn*cosd(phi))*exp(cx*xn+cz*zn);
                end
                Sum_H_mass(p) = Sum_H;
                RCS(p) = 10*log10((2/pi)*Sum_H*conj(Sum_H));

            end

            % график поля в дальней зоне от угла налюденя (не нормированный)
            figure;
            plot(phi_for_graf_DA ,RCS,'r-');
            title('График распределения поля в дальней зоне (MIE MFIE Hpol)'); 
            xlim([0 360])
            xlabel('Угол, град'); 
            ylabel('ЭПР');
            
            % ЗАПИШЕМ
            obj.DA = RCS; 
            obj.phi_for_graf_DA = phi_for_graf_DA;
            obj.Sum_H_mass = Sum_H_mass;
        end
 %%       
        % взаписать токи в файл
        function obj= calculate_DE(obj)
        end
 %%       
        % расчитать невязку
        function obj= write_DA_to_file(obj)
            
            % ВЫТАЩИМ
            phi_for_graf_DA = obj.phi_for_graf_DA;
            DA = obj.DA;
            
            f02 = fopen('DA_MFIE_MIE_Hpol_elips.dat','w');
            for p = 1 : 721
                fprintf(f02, ' %10.5f %10.5f\n', phi_for_graf_DA(p), DA(p));
            end
            fclose(f02);
        end
                %%
        % получить ток -  GETer
        function I = get_I(obj)
            I = obj.I;
        end
        % получить ЭПР -  GETer
        function DA = get_DA(obj)
            DA = obj.DA;
        end
    end
end