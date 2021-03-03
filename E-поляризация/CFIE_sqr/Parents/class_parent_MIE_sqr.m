classdef class_parent_MIE_sqr < handle
%% '''class_parent_MIE_sqr - ������������ ����� ��� MIE sqrt MFIE � EFIE'''

%% ����� ����� (��� ����� ����� ��� ���� �������): 
    % calculate_main_sqr(N)
    % print_sqr()
    % class_parent_MIE_sqr(as, phi_grad)
    % print_sqr()
    % get_I()
    % get_DA()
    % write_I_to_file()
    % write_DA_to_file()
%% ����� �������� �������� (���� ��� ���)
    % ��������� �������: sqr_main (������)
    % ���� � ���: I_real, DA (������)
    % ��� ���: type_MIE (���������)
    % ������������� ������ ������: Sum_E_mass, Zmn, Hi, Ei (������)
    % ��������� �������� ����� � �������: phi_grad, N, as (�����������)
    % ���������: k, gamma (���������)
    % ����������� �� ���� ��������� : delta, dx, n_for_graf, phi_for_graf_DA
    % (������)

%% *** �������� ������ ***
%%  ��� ����������� ��������
    properties 
        sqr_main = struct();  
        I_real = [];
        DA = [];
        % �� ��� ����� ��������� �� ��� 
        Sum_E_mass = [];
        Zmn = [];
        Hi = [];
        Ei = [];
        type_MIE = NaN;
    end
%% ��������� �������� � �������
    properties
        as = NaN;
        phi_grad = NaN;
        N = NaN;
        delta = [];
        dx = [];
        n_for_graf = NaN;
        phi_for_graf_DA = NaN;
    end
    
%%  ���������   
    properties
        k = 2*pi;
        gamma = 1.781072417990;     
    end
    
%% *** ������***
    methods
%% �����������        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % �����������
        function obj = class_parent_MIE_sqr(as, phi_grad)
            obj.as = as;
            obj.phi_grad = phi_grad;
        end
%% ������ ��������        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % ������ ��������� �������
        % ����� ���-�� ����� �� ������ 
        function calculate_main_sqr(obj, N)
            % ������� � �������� N � ������ �� ������������� ����� �����
            obj.N = N;
            
            % ������ ����������� ��������� �� ��-� ������
            % ������� ������
            as = obj.as; %c������� ��������
            k  = obj.k;
            
            % ������ dx
            dx = as/N;
            
            % ��������� �����
            x = zeros(1,N);
            z = zeros(1,N);

            % ��������� �����
            % ������ ��������� ������� �������� ��������� (��� �������)
            % ���������� ��������� ��� ������� �������
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

            % ��������� �� ����������� ��������
            X = [x1 x2 fliplr(x3) x4];
            Z = [z1 z2 z3 fliplr(z4)];

            % ������ ����������� � �������� ���������
            % ����������� * �� �����������
            tx_one = [-1 0 1 0];
            tz_one = [0 -1 0 1];

            %  ����� �����������
            tz = [repmat(tz_one(1),1,N) repmat(tz_one(2),1,N) repmat(tz_one(3),1,N) repmat(tz_one(4),1,N)];
            tx = [repmat(tx_one(1),1,N) repmat(tx_one(2),1,N) repmat(tx_one(3),1,N) repmat(tx_one(4),1,N)];
            
            % ��������� n ��� ������� �������
            %   ���������� ������ ������ � ������� ������� ��������
                a = -90;
                matrix = [cosd(a),-sind(a); sind(a),cosd(a)];

            for i = 1 : 4*N
                 matrix_temp = [tx(i), tz(i)]*matrix;
                 nx(i) = matrix_temp(1);
                 nz(i) = matrix_temp(2);
            end
            
            % ������ ��� ����� ������
            delta = dx/2;
                       
            % ������� ����������� ��������� �� ��-� ������          
            obj.sqr_main.X = X;
            obj.sqr_main.Z = Z;
            obj.sqr_main.tx = tx;
            obj.sqr_main.tz = tz;
            obj.sqr_main.nx = nx;
            obj.sqr_main.nz = nz;
            obj.dx = dx;
            obj.delta = delta;
        end
%% ����� ��������        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ������� ��� ���������� ������� ������� 
        function print_sqr(obj) 
            % �������
            X = obj.sqr_main.X;
            Z = obj.sqr_main.Z;
            X = obj.sqr_main.X;
            Z = obj.sqr_main.Z;
            
            figure;
            plot(X, Z, 'b.');
            title('�������� ������� � ��� ����������� (MAS EFIE Epol)'); 
            legend('class\_MFIE\_Epol\_sqr');
            
        end

%% ������ � ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % �������� ��� -  GETer
        function I = get_I(obj)
            I = obj.I_real;
        end
        % �������� ��� -  GETer
        function DA = get_DA(obj)
            DA = obj.DA;
        end
    end

%%  *** ������***
    methods %(Access = protected)      
%% ������ ������ � ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % �������� ���� � ����
        function obj = write_I_to_file(obj)

            % ������� 
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
        % ��������� �������
        function obj =  write_DA_to_file(obj)
            
            % �������
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

