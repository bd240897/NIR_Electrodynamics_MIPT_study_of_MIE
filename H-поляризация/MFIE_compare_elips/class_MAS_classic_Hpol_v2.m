classdef class_MAS_classic_H_pol < handle
    %%  Уже расчитанные значения
    properties
        I_mas = [];
        I_real = [];
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
    end
    %% Параметры вспомогательного элипса
    properties
        a1 = NaN;
        b2 = NaN;
    end
%%  Константы   
    properties
        k = 2*pi;
    end
    
 %%   
    methods  
%%
        % конструктор
        function obj = class_MAS_classic_Hpol(a,b,phi_grad)
        end
 %%                       
        % расчет вспомогательного эллипсса эллипса
        % выбор кол-ва ВИ
        function obj, ellips_mas = calculate_mas_ellips(N_mas)
        end
%%        
        % расчет основного эллипса
        % выбор кол-ва точек на элипсе 
        function obj, ellips_main  = calculate_main_ellips(N_main)
        end
%%     
        % расчет токов для данного числа источников и точек на элипсе
        function I_real = calculate_I(N_mas, N_main)
            ellips_mas = calculate_mas_ellips(N_mas);
            % достать параметры от сюда
            ellips_main  = calculate_main_ellips(N_main);
            % достать параметры от сюда
        end
%%        
        % записать токи в файл
        function obj = write_I_to_file(obj, N_main)
        end
%%        
        % расчитать ЭПР
        function obj= calculate_DA(obj)
        end
 %%       
        % взаписать токи в файл
        function obj= write_DA_to_file(obj)
        end
 %%       
        % расчитать невязку
        function obj= calculate_DE(obj)
        end
    end
end