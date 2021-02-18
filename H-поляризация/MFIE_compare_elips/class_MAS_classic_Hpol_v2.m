classdef class_MAS_classic_H_pol < handle
    %%  ��� ����������� ��������
    properties
        I_mas = [];
        I_real = [];
        ellips_mas = struct();
        ellips_main = struct();  
    end
%% ��������� ������ � �������
    properties
        a = NaN;
        b = NaN;
        phi_grad = NaN;
        N_mas = NaN;
        N_main = NaN;
    end
    %% ��������� ���������������� ������
    properties
        a1 = NaN;
        b2 = NaN;
    end
%%  ���������   
    properties
        k = 2*pi;
    end
    
 %%   
    methods  
%%
        % �����������
        function obj = class_MAS_classic_Hpol(a,b,phi_grad)
        end
 %%                       
        % ������ ���������������� �������� �������
        % ����� ���-�� ��
        function obj, ellips_mas = calculate_mas_ellips(N_mas)
        end
%%        
        % ������ ��������� �������
        % ����� ���-�� ����� �� ������ 
        function obj, ellips_main  = calculate_main_ellips(N_main)
        end
%%     
        % ������ ����� ��� ������� ����� ���������� � ����� �� ������
        function I_real = calculate_I(N_mas, N_main)
            ellips_mas = calculate_mas_ellips(N_mas);
            % ������� ��������� �� ����
            ellips_main  = calculate_main_ellips(N_main);
            % ������� ��������� �� ����
        end
%%        
        % �������� ���� � ����
        function obj = write_I_to_file(obj, N_main)
        end
%%        
        % ��������� ���
        function obj= calculate_DA(obj)
        end
 %%       
        % ��������� ���� � ����
        function obj= write_DA_to_file(obj)
        end
 %%       
        % ��������� �������
        function obj= calculate_DE(obj)
        end
    end
end