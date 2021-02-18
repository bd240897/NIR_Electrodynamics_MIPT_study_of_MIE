classdef class_MAS_classic_Hpol_v2 < handle
% """ Класс для расчета MAS calssic H-pol """
% Пример использования:
    % a = class_MAS_classic_Hpol_v2(2,1,0);
    % a.calculate_mas_ellips(100);
    % a.calculate_main_ellips(120);
    % a.calculate_I(100,100).write_I_to_file();
    % a.calculate_DA().write_DA_to_file();
    % a.calculate_DE()
    
    %%  Уже расчитанные значения
    properties
        I_mas = [];
        I_real = [];
        DA = [];
        ellips_mas = struct();
        ellips_main = struct();  
    end
%% Параметры элипса и расчета
    properties
        a = NaN;
        b = NaN;
        phi_grad = NaN;
        N_mas = NaN;
        N_main = NaN;
        phi_for_graf_real = NaN;
        phi_for_graf_DA = NaN;
    end
    %% Параметры вспомогательного элипса
    properties
        a1 = NaN;
        b1 = NaN;
    end
%%  Константы   
    properties
        k = 2*pi;
    end
    
 %%   
    methods  
%%
        % конструктор
        function obj = class_MAS_classic_Hpol_v2(a,b,phi_grad)
            obj.a = a;
            obj.b = b;
            obj.phi_grad = phi_grad;
        end
 %%                       
        % расчет вспомогательного эллипсса эллипса
        % выбор кол-ва ВИ
        function obj = calculate_mas_ellips(obj, N_mas)
            % ЗАПИШЕМ в свойства N и бдудем ег оиспользовать везде далее
            obj.N_mas = N_mas;
            
            % вытащи необходимые параметры из СВ-В КЛАССА
            % полуоси элипса
            a = obj.a; %большшая полуось - ОX
            b = obj.b; %меньшая полуось - ОZ
            
            f = sqrt(a^2 - b^2);
            a1 = 1.8;  
            b1 = sqrt(a1^2 - f^2); 
            
            % Генерация точек контура для ВИ
            delta_phi = 360/N_mas;
            for i = 1 : N_mas
                xAS(i) = a1*cosd(delta_phi*(i-1));
                zAS(i) = b1*sind(delta_phi*(i-1));
            end
            
            plot(xAS,zAS, 'xr')
            
            % ЗАПИШЕМ необходимые параметры из СВ-В КЛАССА
            obj.a1 = a1;
            obj.b1 = b1;
            obj.ellips_mas.xAS = xAS;
            obj.ellips_mas.zAS = zAS;
        end
%%        
        % расчет основного эллипса
        % выбор кол-ва точек на элипсе 
        function obj  = calculate_main_ellips(obj, N_main)
            % ЗАПИШЕМ в свойства N и бдудем ег оиспользовать везде далее
            obj.N_main = N_main;
            
            % вытащи необходимые параметры из СВ-В КЛАССА
            % полуоси элипса
            a = obj.a; %большшая полуось - ОX
            b = obj.b; %меньшая полуось - ОZ
            k  = obj.k;

            % Генерация точек элипса
            delta_phi = 360/N_main;
            x = zeros(1,N_main);
            z = zeros(1,N_main);
            
            for i = 1 : N_main
                x(i) = a*cosd(delta_phi*(i-1));
                z(i) = b*sind(delta_phi*(i-1));
            end
            
            plot(x,z)
    
            % Касательная к основному элипсу
            for i = 1 : N_main
                tx(i) = -a^2*z(i)/sqrt(a^4*z(i)^2+b^4*x(i)^2);
                tz(i) =  b^2*x(i)/sqrt(a^4*z(i)^2+b^4*x(i)^2);
            end
            
            % ЗАПИШЕМ необходимые параметры из СВ-В КЛАССА          
            obj.ellips_main.x = x;
            obj.ellips_main.z = z;
            obj.ellips_main.tx = tx;
            obj.ellips_main.tz = tz;
        end
%%     
        % расчет токов для данного числа источников и точек на элипсе
        function obj = calculate_I(obj, N_mas, N_main)
            
            % ВЫТАЩИМ 
            phi_grad = obj.phi_grad;   
            delta_phi = 360/N_main;
            k = obj.k; 
                       
            % генерация точек для ВИ
            obj.calculate_mas_ellips(N_mas);
            % достать параметры от сюда
            xAS = obj.ellips_mas.xAS;
            zAS = obj.ellips_mas.zAS;
           
            % генерация точек освного эллипса
            obj.calculate_main_ellips(N_main);
            % достать параметры от сюда
            x = obj.ellips_main.x;
            z = obj.ellips_main.z;
            tx = obj.ellips_main.tx;
            tz = obj.ellips_main.tz;
            
                        % Изучим один цикл
            % count = 0;
            % пробегаемся по всем точкам оснвого элипса
            for m = 1 : N_main
                xm = x(m);
                zm = z(m);
                txm = tx(m);
                tzm = tz(m);

                %пробегаемся по всем точкам ВИ
                for n = 1 : N_mas
                    xn = xAS(n);
                    zn = zAS(n);

                    % Растояние между m и n
                    dx = xm - xn;
                    dz = zm - zn;
                    r = sqrt(dx^2 + dz^2);

                    Es(m,n) = besselh(1,1,k*r)*(dx*tzm - dz*txm)/r;         

                end
            end

            % Расчет падающего поля
            for m=1:N_main
                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);
                % Определим координаты
                xm = x(m);
                zm = z(m);
                txm = tx(m);
                tzm = tz(m);

            %   Итоговое падающее поел - это суперпозиция полей
                Eix = +1i*sind(phi_grad)*exp(cx*xm+cz*zm);% падающее поле X
                Eiz = -1i*cosd(phi_grad)*exp(cx*xm+cz*zm);% падающее поле Z
                Ei(m,1) = Eix*txm + Eiz*tzm;         
            end


            % Посчитаем токи и построим график
            % расчеток токов ВИ 
            I=Es\Ei;

            for m = 1 : N_main
                xm = x(m);
                zm = z(m);
                txm = tx(m);
                tzm = tz(m);

                %пробегаемся по всем точкам ВИ
                sum_count = 0;
                for n = 1 : N_mas
                    xn = xAS(n);
                    zn = zAS(n);    

                    % Растояние между m и n
                    xmn = (xm - xn);
                    zmn = (zm - zn);
                    r = sqrt(xmn^2 + zmn^2);

                    Hs = I(n)*besselh(0, 1, k*r);

                    sum_count = sum_count + Hs;
                end

                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);
                Hi = exp(cx*xm+cz*zm); 

                I_real(m,1) = sum_count + Hi;
            end

            for i = 1 : N_main
                z_m = z(i);
                x_m = x(i);

                if (x_m <= 0 && z_m < 0) || (x_m >= 0 && z_m < 0)
                    phi_for_graf_real(i) = atan2(z_m, x_m)*180/pi+360;
                else
                    phi_for_graf_real(i) = atan2(z_m, x_m)*180/pi;
                end
            end

            figure;
            fig = plot(phi_for_graf_real, abs(I_real), 'r');
            title('График распределения тока в ВИ'); 
            xlabel('Координаты центра цилиндрика'); 
            xlim([0 360])
            ylabel('Плотность тока');
            
            % ЗАПИШЕМ необходимые параметры из СВ-В КЛАССА
            obj.I_mas = I;
            obj.I_real = I_real; 
            obj.phi_for_graf_real = phi_for_graf_real;
            
        end
%%        
        % записать токи в файл
        function obj = write_I_to_file(obj)
            
            % ВЫТАЩИМ 
            I_real = obj.I_real;
            N_main = obj.N_main;
            phi_for_graf_real = obj.phi_for_graf_real;
            
            name_file = strcat('I_MAS_elips_pres_',  num2str(N_main), '.dat'); % в имея файла добавил число отрезков автоматом 

            f01 = fopen(name_file,'w');
            for n = 1 : N_main
                fprintf(f01,' %10.5f %10.5f\n', phi_for_graf_real(n), abs(I_real(n)));
            end
            fclose(f01);
        end
%%        
        % расчитать ЭПР
        function obj = calculate_DA(obj)
            
            % ВЫТАЩИМ необходимые параметры из СВ-В КЛАССА 
            k = obj.k;
            I = obj.I_mas;
            xAS = obj.ellips_mas.xAS;
            zAS = obj.ellips_mas.zAS;
            N_mas = obj.N_mas;
            
            for p = 1 : 721
                Sum_H = 0;
                phi = (p-1)/2; % угол из формулы для приведения 721 к 360 градусам
                phi_for_graf_DA(p) = (p-1)/2; % угол для графика

                cx = - 1i*k*cosd(phi);
                cz = - 1i*k*sind(phi);

                for n = 1 : N_mas
                   Sum_H = Sum_H + I(n)*exp(cx*xAS(n)+cz*zAS(n)); 
                end

                RCS(p) = 10*log10((2/pi)*Sum_H*conj(Sum_H));

            end

            % график поля в дальней зоне от угла налюденя (не нормированный)
            figure;
            plot(phi_for_graf_DA ,RCS,'r-');
            title('График распределения поля в дальней зоне'); 
            xlim([0 360])
            xlabel('Угол, град'); 
            ylabel('ЭПР');
            
            % ЗАПИШЕМ
            obj.DA = RCS; 
            obj.phi_for_graf_DA = phi_for_graf_DA;
        end
 %%       
        % взаписать токи в файл
        function obj= write_DA_to_file(obj)
            
            % ВЫТАЩИМ 
            DA = obj.DA;
            N_main = obj.N_main;
            phi_for_graf_DA = obj.phi_for_graf_DA;
                     
            name_file = strcat('DA_MAS_elips_pres_',  num2str(N_main), '.dat'); % в имея файла добавил число отрезков автоматом 
            f02 = fopen(name_file,'w');
            for p = 1 : 721
                fprintf(f02,' %10.5f %10.5f\n', phi_for_graf_DA(p), DA(p));
            end
            fclose(f02);
            
        end
%%       
        % расчитать невязку
        function obj = calculate_DE(obj)
            % ВЫТАЩИМ необходимые параметры из СВ-В КЛАССА
            phi_grad = obj.phi_grad;   
            a = obj.a;
            b = obj.b;
            x = obj.ellips_main.x;
            z = obj.ellips_main.z;
            k = obj.k;
            I = obj.I_mas;
            xAS = obj.ellips_mas.xAS;
            zAS = obj.ellips_mas.zAS;
            N_main = obj.N_main;
            N_mas = obj.N_mas;
            
            % Генерация точек элипса
            delta_phi_DE = 360/(2*N_main);
            x = zeros(1,N_main);
            z = zeros(1,N_main);

            for i = 1 : 2*N_main
                x_DE(i) = a*cosd(delta_phi_DE*(i-1));
                z_DE(i) = b*sind(delta_phi_DE*(i-1));
            end

            % plot(x_DE,z_DE)

            % Генерация касательных в новых точках
            for i = 1 : 2*N_main
                tx_DE(i) = -a^2*z_DE(i)/sqrt(a^4*z_DE(i)^2+b^4*x_DE(i)^2);
                tz_DE(i) =  b^2*x_DE(i)/sqrt(a^4*z_DE(i)^2+b^4*x_DE(i)^2);
            end

            for m = 1:2*N_main  
                xm = x_DE(m);
                zm = z_DE(m);
                txm = tx_DE(m);
                tzm = tz_DE(m);

                count(m) = m;
                E_DE_sum = 0;    

                for n = 1:N_mas
                    xn = xAS(n);
                    zn = zAS(n);    

                    % токи ВИ
                    xmn = xm - xn;
                    zmn = zm - zn;
                    r = sqrt(xmn^2 + zmn^2);

                    % рассеянное поле на поверхности в 2*N-1 точках
            %         Hs_DE = I(n)*besselh(0,1,k*r);  
                    % рассеяное касательное электрическое поле в точках для невязка
                    Esx = +1i*I(n)*besselh(1,1,k*r)*zmn/r;
                    Esz = -1i*I(n)*besselh(1,1,k*r)*xmn/r;        
                    Es_DE = Esx*txm + Esz*tzm; 

                    % суммарное поле от ВИ в точке m
                    E_DE_sum = E_DE_sum + Es_DE;

                end

                % коэфиценты
                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);  

                % падаюшее поле на поверхности в 2*N-1 точках
            %     Hi_DE = exp(cx*xm+cz*zm);
                % падающее касательное электрическое поле в точках для невязка
                Eix = +sind(phi_grad)*exp(cx*xm+cz*zm);% падающее поле X
                Eiz = -cosd(phi_grad)*exp(cx*xm+cz*zm);% падающее поле Z
                Ei_DE = Eix*txm + Eiz*tzm;

                % cам расчет невязки
                dE(m,1) =  abs(Ei_DE - E_DE_sum);
            end    

            % график для невязки
            figure;
            plot(count,dE,'r-');
            title('График невязок'); 
            xlabel('Номер точки на поверхности цилиндра'); 
            ylabel('Невязка');
        end
    end
end