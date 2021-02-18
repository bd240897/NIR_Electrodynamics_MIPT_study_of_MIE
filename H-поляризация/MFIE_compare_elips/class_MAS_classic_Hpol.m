classdef class_MAS_classic_Hpol < handle
%%  ��� ����������� ��������
    properties
        I_mas = [];
        I_real = [];
        elips_mas = struct();
        elips_source = struct();
        
    end
%% ��������� ������ � �������
    properties
        a = NaN;
        b = NaN;
        phi_grad = NaN;
        N = NaN;
    end
%%  ���������   
    properties
        k = 2*pi;
    end
%% 
    methods  
        % ����������� - ������� ���������
        function obj = class_MAS_classic_Hpol(a,b,phi_grad)
            obj.a = a;
            obj.b = b;
            obj.phi_grad = phi_grad;
        end
 %%            
        % ��������� �������� ������ � ������� ��� ��
        function obj = generate_ellips(obj, N)
            
            % ������� � �������� N � ������ �� ������������� ����� �����
            obj.N = N;
            
            % ������ ����������� ��������� �� ��-� ������
            % ������� ������
            a = obj.a; %�������� ������� - �X
            b = obj.b; %������� ������� - �Z
            
            f = sqrt(a^2 - b^2);
            
            a1 = 1.8;  
            b1 = sqrt(a1^2 - f^2); 
            
            delta_phi = 360/N;
            k = 2*pi; 
            
            % ��������� ����� ������
            delta_phi = 360/N;
            x = zeros(1,N);
            z = zeros(1,N);
            
            for i = 1 : N
                x(i) = a*cosd(delta_phi*(i-1));
                z(i) = b*sind(delta_phi*(i-1));
            end
            
            plot(x,z)

            % ��������� ����� ������� ��� ��
            for i = 1 : N
                xAS(i) = a1*cosd(delta_phi*(i-1));
                zAS(i) = b1*sind(delta_phi*(i-1));
            end
            
            plot(x,z, xAS,zAS, 'xr')
            
            % ����������� � ��������� ������
            for i = 1 : N
                tx(i) = -a^2*z(i)/sqrt(a^4*z(i)^2+b^4*x(i)^2);
                tz(i) =  b^2*x(i)/sqrt(a^4*z(i)^2+b^4*x(i)^2);
            end
            
            % ������� ����������� ��������� �� ��-� ������
            obj.elips_mas.xAS = xAS;
            obj.elips_mas.zAS = zAS;
            
            obj.elips_source.x = x;
            obj.elips_source.z = z;
            obj.elips_source.tx = tx;
            obj.elips_source.tz = tz;
        end
            
  %%       
        % ������ �����
        function[I_real] = reseach_I(obj, N) 
            
            % ������� ����������� ��������� �� ��-� ������ 
            phi_grad = obj.phi_grad;   
            delta_phi = 360/N;
            k = obj.k  
            x = obj.elips_source.x;
            z = obj.elips_source.z;
            tx = obj.elips_source.tx;
            tz = obj.elips_source.tz;
            xAS = obj.elips_mas.xAS;
            zAS = obj.elips_mas.zAS;
            
            
            % ������ ���� ����
            % count = 0;
            % ����������� �� ���� ������ ������� ������
            for m = 1 : N
                xm = x(m);
                zm = z(m);
                txm = tx(m);
                tzm = tz(m);

                %����������� �� ���� ������ ��
                for n = 1 : N
                    xn = xAS(n);
                    zn = zAS(n);

                    % ��������� ����� m � n
                    dx = xm - xn;
                    dz = zm - zn;
                    r = sqrt(dx^2 + dz^2);

                    Es(m,n) = besselh(1,1,k*r)*(dx*tzm - dz*txm)/r;         

                end
            end

            % ������ ��������� ����
            for m=1:N
                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);
                % ��������� ����������
                xm = x(m);
                zm = z(m);
                txm = tx(m);
                tzm = tz(m);

            %   �������� �������� ���� - ��� ������������ �����
                Eix = +1i*sind(phi_grad)*exp(cx*xm+cz*zm);% �������� ���� X
                Eiz = -1i*cosd(phi_grad)*exp(cx*xm+cz*zm);% �������� ���� Z
                Ei(m,1) = Eix*txm + Eiz*tzm;         
            end


            % ��������� ���� � �������� ������
            % �������� ����� �� 
            I=Es\Ei;

            for m = 1 : N
                xm = x(m);
                zm = z(m);
                txm = tx(m);
                tzm = tz(m);

                %����������� �� ���� ������ ��
                sum_count = 0;
                for n = 1 : N
                    xn = xAS(n);
                    zn = zAS(n);    

                    % ��������� ����� m � n
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

            % ������� ������� ��� ������� � ��������
            % ����� ���� ����� � �������� �������
            begin_phi = delta_phi/2;
            phi_for_graf = [];
            for i = 1 : N 
                count = delta_phi*(i-1);
                phi_for_graf(i) = count;
            end

            figure;
            fig = plot(phi_for_graf, abs(I_real), 'r');
            title('������ ������������� ���� � ��'); 
            xlabel('���������� ������ ����������'); 
            xlim([0 360])
            ylabel('��������� ����');
            
            % ������� ����������� ��������� �� ��-� ������
            obj.I_mas = I;
            obj.I_real = I_real;           
        end    
        
 %%        
        % ������ ���
        function[] = reseach_DA(obj)
            
            % ������� ����������� ��������� �� ��-� ������ 
            k = obj.k;
            I = obj.I_mas;
            xAS = obj.elips_mas.xAS;
            zAS = obj.elips_mas.zAS;
            N = obj.N;
            
            for p = 1 : 721
                Sum_H = 0;
                phi = (p-1)/2; % ���� �� ������� ��� ���������� 721 � 360 ��������
                phi_for_graf_DA(p) = (p-1)/2; % ���� ��� �������

                cx = - 1i*k*cosd(phi);
                cz = - 1i*k*sind(phi);

                for n = 1 : N
                   Sum_H = Sum_H + I(n)*exp(cx*xAS(n)+cz*zAS(n)); 
                end

                RCS(p) = 10*log10((2/pi)*Sum_H*conj(Sum_H));

            end

            % ������ ���� � ������� ���� �� ���� �������� (�� �������������)
            figure;
            plot(phi_for_graf_DA ,RCS,'r-');
            title('������ ������������� ���� � ������� ����'); 
            xlim([0 360])
            xlabel('����, ����'); 
            ylabel('���');
        end
  %%        
        % ������ �������
        function[] = reseach_DE(obj)
            
            % ������� ����������� ��������� �� ��-� ������
            phi_grad = obj.phi_grad;   
            N = obj.N;
            delta_phi = 360/N;
            x = obj.elips_source.x;
            z = obj.elips_source.z;
            k = obj.k;
            I = obj.I_mas;
            xAS = obj.elips_mas.xAS;
            zAS = obj.elips_mas.zAS;
            a = obj.a;
            b = obj.b;
            
            % ��������� ����� ������
            delta_phi_DE = 360/(2*N);
            x = zeros(1,N);
            z = zeros(1,N);

            for i = 1 : 2*N
                x_DE(i) = a*cosd(delta_phi_DE*(i-1));
                z_DE(i) = b*sind(delta_phi_DE*(i-1));
            end

            % plot(x_DE,z_DE)

            % ��������� ����������� � ����� ������
            for i = 1 : 2*N
                tx_DE(i) = -a^2*z_DE(i)/sqrt(a^4*z_DE(i)^2+b^4*x_DE(i)^2);
                tz_DE(i) =  b^2*x_DE(i)/sqrt(a^4*z_DE(i)^2+b^4*x_DE(i)^2);
            end

            for m = 1:2*N  
                xm = x_DE(m);
                zm = z_DE(m);
                txm = tx_DE(m);
                tzm = tz_DE(m);

                count(m) = m;
                E_DE_sum = 0;    

                for n = 1:N
                    xn = xAS(n);
                    zn = zAS(n);    

                    % ���� ��
                    xmn = xm - xn;
                    zmn = zm - zn;
                    r = sqrt(xmn^2 + zmn^2);

                    % ���������� ���� �� ����������� � 2*N-1 ������
            %         Hs_DE = I(n)*besselh(0,1,k*r);  
                    % ��������� ����������� ������������� ���� � ������ ��� �������
                    Esx = +1i*I(n)*besselh(1,1,k*r)*zmn/r;
                    Esz = -1i*I(n)*besselh(1,1,k*r)*xmn/r;        
                    Es_DE = Esx*txm + Esz*tzm; 

                    % ��������� ���� �� �� � ����� m
                    E_DE_sum = E_DE_sum + Es_DE;

                end

                % ����������
                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);  

                % �������� ���� �� ����������� � 2*N-1 ������
            %     Hi_DE = exp(cx*xm+cz*zm);
                % �������� ����������� ������������� ���� � ������ ��� �������
                Eix = +sind(phi_grad)*exp(cx*xm+cz*zm);% �������� ���� X
                Eiz = -cosd(phi_grad)*exp(cx*xm+cz*zm);% �������� ���� Z
                Ei_DE = Eix*txm + Eiz*tzm;

                % c�� ������ �������
                dE(m,1) =  abs(Ei_DE - E_DE_sum);
            end    

            % ������ ��� �������
            figure;
            plot(count,dE,'r-');
            title('������ �������'); 
            xlabel('����� ����� �� ����������� ��������'); 
            ylabel('�������');
        end    
        
    end    
end