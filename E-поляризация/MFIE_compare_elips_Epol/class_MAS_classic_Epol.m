classdef class_MAS_classic_Epol < handle
%% """ ����� ��� ������� MAS classic E-pol """
% 19-02-21

% ������ �������������:
    % a = 2;
    % b = 1;
    % phi_grad = 60;
    % N_mas = 400;
    % N_main= 100;
    % 
    % calc_1 = class_MAS_classic_Epol(a,b,phi_grad);
    % calc_1.calculate_I(N_mas,N_main).write_I_to_file();
    % calc_1.calculate_DA().write_DA_to_file();
    % calc_1.calculate_DE();
    % I = calc_1.get_I();
    % DA = calc_1.get_DA();
    
    %%  ��� ����������� ��������
    properties
        ellips_mas = struct();
        ellips_main = struct();  
        I_mas = [];
        I_real = [];
        Sum_E_mass = [];
        DA = [];
        dE = []
    end
%% ��������� ������ � �������
    properties
        a = NaN;
        b = NaN;
        phi_grad = NaN;
        N_mas = NaN;
        N_main = NaN;
        phi_for_graf_real = NaN;
        phi_for_graf_DA = NaN;
    end
    %% ��������� ���������������� ������
    properties
        a1 = NaN;
        b1 = NaN;
    end
%%  ���������   
    properties
        k = 2*pi;
    end
    
 %%   
    methods  
%%
        % �����������
        function obj = class_MAS_classic_Epol(a,b,phi_grad)
            obj.a = a;
            obj.b = b;
            obj.phi_grad = phi_grad;
        end
 %%                       
        % ������ ���������������� �������� �������
        % ����� ���-�� ��
        function obj = calculate_mas_ellips(obj, N_mas)
            % ������� � �������� N � ������ �� ������������� ����� �����
            obj.N_mas = N_mas;
            
            % ������ ����������� ��������� �� ��-� ������
            % ������� ������
            a = obj.a; %�������� ������� - �X
            b = obj.b; %������� ������� - �Z
            
            f = sqrt(a^2 - b^2);
            a1 = 1.8;  
            b1 = sqrt(a1^2 - f^2); 
            
            % ��������� ����� ������� ��� ��
            delta_phi = 360/N_mas;
            for i = 1 : N_mas
                xAS(i) = a1*cosd(delta_phi*(i-1));
                zAS(i) = b1*sind(delta_phi*(i-1));
            end
            
%             plot(xAS,zAS, 'xr')
%             title('(MAS classic Epol)');
            
            % ������� ����������� ��������� �� ��-� ������
            obj.a1 = a1;
            obj.b1 = b1;
            obj.ellips_mas.xAS = xAS;
            obj.ellips_mas.zAS = zAS;
        end
%%        
        % ������ ��������� �������
        % ����� ���-�� ����� �� ������ 
        function obj  = calculate_main_ellips(obj, N_main)
            % ������� � �������� N � ������ �� ������������� ����� �����
            obj.N_main = N_main;
            
            % ������ ����������� ��������� �� ��-� ������
            % ������� ������
            a = obj.a; %�������� ������� - �X
            b = obj.b; %������� ������� - �Z
            k  = obj.k;

            % ��������� ����� ������
            delta_phi = 360/N_main;
            x = zeros(1,N_main);
            z = zeros(1,N_main);
            
            for i = 1 : N_main
                x(i) = a*cosd(delta_phi*(i-1));
                z(i) = b*sind(delta_phi*(i-1));
            end
            
%             plot(x,z)
%             title('(MAS classic Epol)');
%             
            % ����������� � ��������� ������
            for i = 1 : N_main
                tx(i) = -a^2*z(i)/sqrt(a^4*z(i)^2+b^4*x(i)^2);
                tz(i) =  b^2*x(i)/sqrt(a^4*z(i)^2+b^4*x(i)^2);
            end
            
            % ������� ����������� ��������� �� ��-� ������          
            obj.ellips_main.x = x;
            obj.ellips_main.z = z;
            obj.ellips_main.tx = tx;
            obj.ellips_main.tz = tz;
        end
%% 
        % ������� ��� ���������� ������� �������� 
        function print_ellips(obj) 
            % �������
            xAS = obj.ellips_mas.xAS;
            zAS = obj.ellips_mas.zAS;
            x = obj.ellips_main.x;
            z = obj.ellips_main.z;
            
            figure;
            plot(x,z, xAS,zAS, 'xr');
            title('�������� � ��������������� ������ (MAS Epol)'); 

        end

%%     
        % ������ ����� ��� ������� ����� ���������� � ����� �� ������
        function obj = calculate_I(obj, N_mas, N_main)
            
            % ������� 
            phi_grad = obj.phi_grad;   
            delta_phi = 360/N_main;
            k = obj.k; 
                       
            % ��������� ����� ��� ��
            obj.calculate_mas_ellips(N_mas);
            % ������� ��������� �� ����
            xAS = obj.ellips_mas.xAS;
            zAS = obj.ellips_mas.zAS;
           
            % ��������� ����� ������� �������
            obj.calculate_main_ellips(N_main);
            % ������� ��������� �� ����
            x = obj.ellips_main.x;
            z = obj.ellips_main.z;
            tx = obj.ellips_main.tx;
            tz = obj.ellips_main.tz;
            
            % �������� ������
            obj.print_ellips()
            
            % ������ ���� ����
            % count = 0;
            % ����������� �� ���� ������ ������� ������
            for m = 1 : N_main
                xm = x(m);
                zm = z(m);

                %����������� �� ���� ������ ��
                for n = 1 : N_mas
                    xn = xAS(n);
                    zn = zAS(n);

                    % ��������� ����� m � n
                    xmn = xm - xn;
                    zmn = zm - zn;
                    r = sqrt(xmn^2 + zmn^2);

                    % ���������� ����
                    Es(m,n) = besselh(0,1,k*r);         

                end
            end


            % ������ ��������� ����
            for m=1:N_main
                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);
                % ��������� ����������
                xm = x(m);
                zm = z(m);

                % ���������� ����
                Ei(m,1) = exp(cx*xm+cz*zm);         
            end


            % ��������� ���� �� (�� �������� ����!)
            % �������� ����� �� 
            I=Es\Ei;

            % ������ �������� ����

            for m = 1 : N_main
                xm = x(m);
                zm = z(m);
                txm = tx(m);
                tzm = tz(m);

                %����������� �� ���� ������ ��
                Hs_sum = 0;
                for n = 1 : N_mas
                    xn = xAS(n);
                    zn = zAS(n);    

                    % ��������� ����� m � n
                    xmn = (xm - xn);
                    zmn = (zm - zn);
                    r = sqrt(xmn^2 + zmn^2);

                    % ����������� ���������� ����
                    Hsx = +1i*I(n)*besselh(1,1,k*r)*zmn/r;
                    Hsz = -1i*I(n)*besselh(1,1,k*r)*xmn/r;

                    Hs_n = Hsx*txm + Hsz*tzm;  

                    Hs_sum = Hs_sum + Hs_n;
                end

                % ����������� �������� ����
                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);

                Hix = +sind(phi_grad)*exp(cx*xm+cz*zm);% �������� ���� X
                Hiz = -cosd(phi_grad)*exp(cx*xm+cz*zm);% �������� ���� Z
                Hi = Hix*txm + Hiz*tzm;

                I_real(m,1) = Hs_sum + Hi;
            end

            % ����� ���� ��0 ������� �� �������� � �������� - ���� �������� ��������� ���� � ��
%             I_real = I_real*(1/(120*pi));

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
            title('������ ������������� ���� � �� (MAS classic Epol)'); 
            xlabel('����'); 
            xlim([0 360])
            ylabel('��������� ����');
            
            % ������� ����������� ��������� �� ��-� ������
            obj.I_mas = I;
            obj.I_real = I_real; 
            obj.phi_for_graf_real = phi_for_graf_real;
            
        end
%%        
        % �������� ���� � ����
        function obj = write_I_to_file(obj)
            
            % ������� 
            I_real = obj.I_real;
            N_main = obj.N_main;
            phi_for_graf_real = obj.phi_for_graf_real;
            
            name_file = strcat('I_MAS_classic_Epol_elips_',  num2str(N_main), '.dat'); % � ���� ����� ������� ����� �������� ��������� 

            f01 = fopen(name_file,'w');
            for n = 1 : N_main
                fprintf(f01,' %10.5f %10.5f\n', phi_for_graf_real(n), abs(I_real(n)));
            end
            fclose(f01);
        end
%%        
        % ��������� ���
        function obj = calculate_DA(obj)
            
            % ������� ����������� ��������� �� ��-� ������ 
            k = obj.k;
            I = obj.I_mas;
            xAS = obj.ellips_mas.xAS;
            zAS = obj.ellips_mas.zAS;
            N_mas = obj.N_mas;
            
            for p = 1 : 721

                Sum_E = 0;
                phi = (p-1)/2; % ���� �� ������� ��� ���������� 721 � 360 ��������
                phi_for_graf_DA(p) = (p-1)/2; % ���� ��� �������

                cx = - 1i*k*cosd(phi);
                cz = - 1i*k*sind(phi);

                for n = 1 : N_mas

                   % ���������� ����������
                   xn = xAS(n);
                   zn = zAS(n);

                   % ��� �������� �� ���� ��!!!
                   Sum_E = Sum_E + I(n)*exp(cx*xn+cz*zn); 
                end
                
                Sum_E_mass(p) = Sum_E;
                RCS(p) = 10*log10((2/pi)*Sum_E*conj(Sum_E));

            end

            % ������ ���� � ������� ���� �� ���� �������� (�� �������������)
            figure;
            plot(phi_for_graf_DA ,RCS,'r-');
            title('������ ������������� ���� � ������� ���� (MAS classic Epol)'); 
            xlim([0 360])
            xlabel('����, ����'); 
            ylabel('���');
            
            % �������
            obj.DA = RCS; 
            obj.phi_for_graf_DA = phi_for_graf_DA;
            obj.Sum_E_mass = Sum_E_mass;
        end
 %%       
        % ��������� ���� � ����
        function obj= write_DA_to_file(obj)
            
            % ������� 
            DA = obj.DA;
            N_main = obj.N_main;
            phi_for_graf_DA = obj.phi_for_graf_DA;
                     
            name_file = strcat('DA_MAS_classic_Epol_elips_',  num2str(N_main), '.dat'); % � ���� ����� ������� ����� �������� ��������� 
            f02 = fopen(name_file,'w');
            for p = 1 : 721
                fprintf(f02,' %10.5f %10.5f\n', phi_for_graf_DA(p), DA(p));
            end
            fclose(f02);
            
        end
%%       
        % ��������� �������
        function obj = calculate_DE(obj)
            % ������� ����������� ��������� �� ��-� ������
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
            
            % ��������� ����� ������
            delta_phi_DE = 360/(2*N_main);
            x = zeros(1,N_main);
            z = zeros(1,N_main);

            for i = 1 : 2*N_main
                x_DE(i) = a*cosd(delta_phi_DE*(i-1));
                z_DE(i) = b*sind(delta_phi_DE*(i-1));
            end

            % plot(x_DE,z_DE)

            % ��������� ����������� � ����� ������
            for i = 1 : 2*N_main
                tx_DE(i) = -a^2*z_DE(i)/sqrt(a^4*z_DE(i)^2+b^4*x_DE(i)^2);
                tz_DE(i) =  b^2*x_DE(i)/sqrt(a^4*z_DE(i)^2+b^4*x_DE(i)^2);
            end

            for m = 1:2*N_main   
                xm = x_DE(m);
                zm = z_DE(m);

                count(m) = m;
                E_DE_sum = 0;    

                for n = 1 : N_mas
                    xn = xAS(n);
                    zn = zAS(n);    

                    % ���� ��
                    xmn = xm - xn;
                    zmn = zm - zn;
                    r = sqrt(xmn^2 + zmn^2);

                    % ���������� ���� �� ����������� � 2*N-1 ������
                    Es_DE = I(n)*besselh(0,1,k*r);  

                    % ��������� ���� �� �� � ����� m
                    E_DE_sum = E_DE_sum + Es_DE;

                end

                % ����������
                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);  

                % �������� ���� �� ����������� � 2*N-1 ������
                Ei_DE = exp(cx*xm+cz*zm);

                % c�� ������ �������
                dE(m,1) =  abs(Ei_DE - E_DE_sum);
            end     

            % ������ ��� �������
            figure;
            plot(count,dE,'r-');
            title('������ ������� (MAS classic Epol)'); 
            xlabel('����� ����� �� ����������� ��������'); 
            ylabel('�������');
            
            % ������� 
            obj.dE = dE;
        end
        
        %%
        % �������� ��� -  GETer
        function I = get_I(obj)
            I = obj.I_real;
        end
        % �������� ��� -  GETer
        function DA = get_DA(obj)
            DA = obj.DA;
        end
        
    end
end