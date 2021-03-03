classdef class_parent_MIE_sqr < handle
%% '''class_parent_MIE_sqr - Родительский класс для MIE sqrt MFIE и EFIE'''

%% Класс умеет (Это общие куски для всех классов): 
    % calculate_main_sqr(N)
    % print_sqr()
    % class_parent_MIE_sqr(as, phi_grad)
    % print_sqr()
    % get_I()
    % get_DA()
    % write_I_to_file()
    % write_DA_to_file()
%% Класс содержит свойства (поля дли них)
    % параметры эллипса: sqr_main (расчет)
    % Токи и ЭПР: I_real, DA (расчет)
    % Тип МИУ: type_MIE (заполнены)
    % Промежуточные данные вывода: Sum_E_mass, Zmn, Hi, Ei (расчет)
    % Параметры падающей волны и расчета: phi_grad, N, as (конструктор)
    % Константы: k, gamma (заполнены)
    % Расчитанные по ходу перменные : delta, dx, n_for_graf, phi_for_graf_DA
    % (расчет)

%% *** СВОЙСТВА КЛАССА ***
%%  Уже расчитанные значения
    properties 
        sqr_main = struct();  
        I_real = [];
        DA = [];
        % то что будем дсотавать из вне 
        Sum_E_mass = [];
        Zmn = [];
        Hi = [];
        Ei = [];
        type_MIE = NaN;
    end
%% Параметры квадрата и расчета
    properties
        as = NaN;
        phi_grad = NaN;
        N = NaN;
        delta = [];
        dx = [];
        n_for_graf = NaN;
        phi_for_graf_DA = NaN;
    end
    
%%  Константы   
    properties
        k = 2*pi;
        gamma = 1.781072417990;     
    end
    
%% *** МЕТОДЫ***
    methods
%% КОНСТРУКТОР        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % конструктор
        function obj = class_parent_MIE_sqr(as, phi_grad)
            obj.as = as;
            obj.phi_grad = phi_grad;
        end
%% РАСЧЕТ КВАДРАТА        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % расчет основного эллипса
        % выбор кол-ва точек на элипсе 
        function calculate_main_sqr(obj, N)
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

            % Расчет касательных к отрезкам разибения
            % касательные * до копирования
            tx_one = [-1 0 1 0];
            tz_one = [0 -1 0 1];

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
%% ВЫВОД КВАДРАТА        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            legend('class\_MFIE\_Epol\_sqr');
            
        end

%% ГЕТЕРЫ И СЕТЕОР
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % получить ток -  GETer
        function I = get_I(obj)
            I = obj.I_real;
        end
        % получить ЭПР -  GETer
        function DA = get_DA(obj)
            DA = obj.DA;
        end
    end

%%  *** МЕТОДЫ***
    methods %(Access = protected)      
%% МЕТОДЫ ЗАПИСИ В ФАЙЛ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % записать токи в файл
        function obj = write_I_to_file(obj)

            % ВЫТАЩИМ 
            I_real = obj.I_real;
            N = obj.N;
            n_for_graf = obj.n_for_graf;
            dx = obj.dx;

            name_file = strcat('I_', obj.type_MIE ,'_MIE_Epol_sqrt_', num2str(N),'.dat');

            f01 = fopen(name_file,'w');
            for n=1:N
                fprintf(f01,' %10.5f %10.5f\n', n_for_graf(n)*dx, abs(I_real(n)));
            end         
            fclose(f01);         
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % расчитать невязку
        function obj =  write_DA_to_file(obj)
            
            % ВЫТАЩИМ
            phi_for_graf_DA = obj.phi_for_graf_DA;
            DA = obj.DA;
            N = obj.N;
            
            name_file = strcat('DA_', obj.type_MIE ,'_MIE_Epol_sqrt_', num2str(N),'.dat');

            f02 = fopen(name_file,'w');
            for p = 1 : 721
                fprintf(f02, ' %10.5f %10.5f\n', phi_for_graf_DA(p), DA(p));
            end
            fclose(f02);
        end
    end
end

