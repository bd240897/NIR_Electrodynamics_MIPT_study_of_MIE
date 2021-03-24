classdef class_MIE_diel_circl_E_pol < handle
%% """����� ��� MIE EFIE E-pol sqr"""
% 02-03-21

% ������ �������������:
    % as = 1;
    % phi_grad = 60;
    % N= 100;
    %     
    % calc_1 = class_MIE_diel_circl_E_pol(e1, mu1, phi_i_grad);
    % calc_1.calculate_main_circl(a, N)
    % calc_1.calculate_I(a, N);
    % calc_1.calculate_DA();
    
    %%  ��� ����������� ��������
    properties 
        circl_main = struct();  
        I_real = [];
        DA = [];
        % �� ��� ����� ��������� �� ��� 
        Sum_E_mass = [];
        Z = [];
        E = [];
        j1 = [];
        j2 = [];
    end
%% ��������� ������ � �������
    properties
        %  1 - ������� ������ 
        k0 = 2*pi;
        e1 = NaN;
        mu1 = NaN;
        k1 = NaN;
        eta1 = NaN;      
        % ��������� �����-2
        k2 = NaN;
        eta2 = NaN;;

        % ��������� ���������
        a = NaN;  % ������ �����
        N = NaN; % ���������� ����� �� �����
        phi_i_grad = NaN; % ���� ������� [����]
        phi_i_rad = NaN; % ���� ������� [���]
        dx = NaN;
        phi_circl_rad = NaN;
        phi_circl_grad = NaN;
        phi_for_graf_DA = NaN;
    end
    %% ��������� ���������������� ������
    properties
        
    end
%%  ���������   
    properties
        gamma = 1.781072417990;     
    end
    
 %%   
    methods  
%%
        % �����������
        function obj = class_MIE_diel_circl_E_pol(e1, mu1, phi_i_grad)
            % �������
            k0 = obj.k0;
            obj.k2 = k0;
            obj.eta2 = 120*pi;
            
            obj.e1 = e1;
            obj.mu1 = mu1;
            obj.k1 = k0*sqrt(mu1*e1);
            obj.eta1 = 120*pi*sqrt(mu1/e1);
            
            obj.phi_i_grad = phi_i_grad;
            obj.phi_i_rad = 2*pi/360*phi_i_grad;
            
            
        end
%%        
        % ������ ��������� �������
        % ����� ���-�� ����� �� ������ 
        function calculate_main_circl(obj, a, N)
            % ������� � �������� N � ������ �� ������������� ����� �����
            obj.N = N;
            
            % ������ ����������� ��������� �� ��-� ������
            
            %  ���������� �����
            dphi_grad = 360/N; % ������� ��� �� ����� [����]
            dphi_rad =  2*pi/360*dphi_grad; % ������� ��� �� ����� [���]
            dx = 2*a*sin(dphi_rad/2);
            dx2 = dx/2; 
            % ����� � �������� ������� ��������� (����� �����)
            % ��������� �������� ����� � ������� ��������
            phi_circl_rad = zeros(1, N);
            phi_circl_grad = zeros(1, N);
            a_for_graph = zeros(1, N);
            for i = 1 : N
                phi_circl_rad(i) = (i-1)*dphi_rad;
                phi_circl_grad(i) = (i-1)*dphi_grad;
                a_for_graph(i) = a;
            end
%             ��������� �����
            dphi_grad = 360/N;
            x = zeros(1,N);
            z = zeros(1,N);

            for i = 1 : N
                z(i) = a*sind(dphi_grad*(i-1));
                x(i) = a*cosd(dphi_grad*(i-1));
            end
            
            figure;
            plot(x,z)

%             ��������� t ��� ������� �������
            for i = 1 : N
                if i == N
                    tx(i) = -((x(1)-x(i)))/dx;
                    tz(i) = -((z(1)-z(i)))/dx;
                else
                    tx(i) = -(x(i+1)-x(i))/dx;
                    tz(i) = -(z(i+1)-z(i))/dx;
                end
            end
            
%             ��������� n ��� ������� �������
            %   ���������� ������ ������ � ������� ������� ��������
                alpha = -90;
                matrix = [cosd(alpha),-sind(alpha); sind(alpha),cosd(alpha)];

            for i = 1 : N
                 matrix_temp = [tx(i), tz(i)]*matrix;
                 nx(i) = matrix_temp(1);
                 nz(i) = matrix_temp(2);
            end

            % ���������� ������� ��������
            for i = 1 : N 
                if i == N
                    x_midle(i) = ((x(1)+x(i)))/2;
                    z_midle(i) = ((z(1)+z(i)))/2;
                else
                    x_midle(i) = (x(i+1)+x(i))/2;
                    z_midle(i) = (z(i+1)+z(i))/2;
                end    
            end
            
            figure;
            plot(x,z,x_midle, z_midle, 'rx');

            % ������� ����������� ��������� �� ��-� ������          
            obj.circl_main.x_midle = x_midle;
            obj.circl_main.z_midle = z_midle;
            obj.circl_main.tx = tx;
            obj.circl_main.tz = tz;
            obj.circl_main.nx = nx;
            obj.circl_main.nz = nz;
            obj.dx = dx;
            obj.phi_circl_rad = phi_circl_rad;
            obj.phi_circl_grad = phi_circl_grad;
            obj.a = a;
        end
%% 
%         % ������� ��� ���������� ������� ������� 
%         function print_sqr(obj) 
%             % �������
%             X = obj.sqr_main.X;
%             Z = obj.sqr_main.Z;
%             X = obj.sqr_main.X;
%             Z = obj.sqr_main.Z;
%             
%             figure;
%             plot(X, Z, 'b.');
%             title('�������� ������� � ��� ����������� (MAS EFIE Epol)'); 
%             legend('class\_MFIE\_Epol\_sqr');
%             
%         end
%%     
        % ������ ����� ��� ������� ����� ���������� � ����� �� ������
        function obj = calculate_I(obj, a, N)
            
            % �������
%             N = obj.N;
%             a = obj.a; 
            
            % ��������� ����� ��� �������
            obj.calculate_main_circl(a, N);
            
            % ������� ��������� �� ����
            x_midle = obj.circl_main.x_midle;
            z_midle = obj.circl_main.z_midle;
            tx = obj.circl_main.tx;
            tz = obj.circl_main.tz;
            nx = obj.circl_main.nx; 
            nz = obj.circl_main.nz; 
            dx = obj.dx;
            phi_i_grad = obj.phi_i_grad;
            phi_i_rad = obj.phi_i_rad;
            k1  = obj.k1;
            k2 = obj.k2;
            eta1 = obj.eta1;
            eta2 = obj.eta2;
            phi_circl_rad = obj.phi_circl_rad;
            phi_circl_grad = obj.phi_circl_grad;
            gamma = obj.gamma;
            
            
            
            % �������� ������
%             obj.print_sqr()
                               
            % ������ ���� ����
            % 3 - c�� ���� ������� 
            % ������ ��� �������� ��������� ������� 
            Z = zeros(2*N,2*N);
            Y = zeros(2*N,2*N);

            for m  = 1 : N
                xm = x_midle(m);
                zm = z_midle(m);
                txm = tx(m);
                tzm = tz(m);
                nxm = nx(m);
                nzm = nz(m);

                for n = 1 : N
                    xn = x_midle(n);
                    zn = z_midle(n);
                    txn = tx(n);
                    tzn = tz(n);
                    nxn = nx(n);
                    nzn = nz(n);

                    % ��������� ����� ������� ��������� �2 � ������ ���������� (� ���)        
                    x_mO2 = xm - xn;
                    z_mO2 = zm - zn;

                    % ������� ���� ������ �� ����� ��
                    % ���������� ���� ����� (� ���)
                    umn = x_mO2*txn + z_mO2*tzn;
                    wmn = x_mO2*nxn + z_mO2*nzn;


                    % ��������� ������ 
                    xmn = xm - xn;
                    zmn = zm - zn;
                    r = sqrt(xmn^2 + zmn^2);

                    umn_plus = umn + dx/2;
                    umn_minus = umn - dx/2; 
                    r_minus = sqrt(umn_minus^2 + wmn^2);
                    r_plus = sqrt(umn_plus^2 + wmn^2);
                    rmn = sqrt(umn^2 + wmn^2);

                    if n == m
                        % ������������� ����
                        E1y = -k1*eta1*dx*1/4*(1 + 1i*2/pi*(log(gamma*k1*dx/4)-1));
                        E2y = -k2*eta2*dx*1/4*(1 + 1i*2/pi*(log(gamma*k2*dx/4)-1));

                        % ��������� ����
                        H1u = -1/2;
                        H2u = +1/2;
                        H1w = 0;
                        H2w = 0;

                    elseif abs(n-m) <= -1.5

                        % ������������� ����
                        E1y = k1*dx + k1*2*1i/pi * (wmn*(atan(umn_plus/wmn) - atan(umn_minus/wmn))...
                        + umn_plus*log(gamma*k1/2*sqrt(umn_plus^2+wmn^2))...
                        - umn_minus*log(gamma*k1/2*sqrt(umn_minus^2+wmn^2)) - dx);
                        E1y = -E1y * eta1 * 1/4;

                        E2y = k2*dx + k2*2*1i/pi * (wmn*(atan(umn_plus/wmn) - atan(umn_minus/wmn))...
                        + umn_plus*log(gamma*k2/2*sqrt(umn_plus^2+wmn^2))...
                        - umn_minus*log(gamma*k2/2*sqrt(umn_minus^2+wmn^2)) - dx);
                        E2y = -E2y * eta2 * 1/4;

                        % ��������� ���� 
                        H1u = 1i/8*wmn*k1^2*dx + 1/(2*pi)*(atan(umn_plus/wmn) - atan(umn_minus/wmn)); % 
                        H2u = 1i/8*wmn*k2^2*dx + 1/(2*pi)*(atan(umn_plus/wmn) - atan(umn_minus/wmn)); % delta?

                        H1w = 1i/4*(obj.H(0,k1*r_plus) - obj.H(0,k1*r_minus));
                        H2w = 1i/4*(obj.H(0,k2*r_plus) - obj.H(0,k2*r_minus));
                    else
                        % ������������� ����
                        E1y = -k1*eta1*1/4*dx*obj.H(0,k1*rmn); % delta?
                        E2y = -k2*eta2*1/4*dx*obj.H(0,k2*rmn); % delta?

                        % ��������� ����     
                        H1u = 1i/4*k1*dx*(wmn/rmn*obj.H(1,k1*rmn)); % delta?
                        H2u = 1i/4*k2*dx*(wmn/rmn*obj.H(1,k2*rmn)); % delta?

                        H1w = 1i/4*(obj.H(0,k1*r_plus) - obj.H(0,k1*r_minus));
                        H2w = 1i/4*(obj.H(0,k2*r_plus) - obj.H(0,k2*r_minus));
                    end
                    % �������� �������       
                    E1mn = E1y;
                    E2mn = -E2y;

                    dot_tm_tn = (txm*txn+tzm*tzn);
                    dot_tm_nn = (txm*nxn+tzm*nzn);

                    H1mn = dot_tm_tn*H1u + dot_tm_nn*H1w;
                    H2mn = -(dot_tm_tn*H2u + dot_tm_nn*H2w);

                    % ������ � ������ ��� ���������� ������� 
                    Z(m,n) = E1mn;
                    Z(m,n+N) = E2mn;
                    Z(m+N,n) = H1mn;
                    Z(m+N,n+N) = H2mn;
                end
            end

            E = zeros(2*N,1);
            % �������� ���� 
            for m = 1 : N
                xm = x_midle(m);
                zm = z_midle(m);
                txm = tx(m);
                tzm = tz(m);

                % c����������� �����
                cx = -1i*k2*cosd(phi_i_grad);
                cz = -1i*k2*sind(phi_i_grad);

                Eiy = exp(cx*xm+cz*zm);
                Hix = 1/eta2*sind(phi_i_grad) * Eiy;
                Hiz = -1/eta2*cosd(phi_i_grad) * Eiy;

                % �������� �������  
                Ei = Eiy;
                Hi = txm*Hix + tzm*Hiz;

                % ������ � ������ ��� ���������� ������� 
                E(m) = Ei;
                E(m+N) = Hi;
            end    

            % ��������� ����
            % ��� ������ ����� 2 ����
            I = Z\E;

            % ������� ���� ��� �� �������
            j1 = I(1:N);
            j2 = I(N+1:2*N);

            % �������� ��� ������� 
            figure;
            plot(phi_circl_grad,  abs(j1), 'r',phi_circl_grad, abs(j2), 'b');
            figure;
            % �������� ������ �����
            polarplot(phi_circl_rad, abs(j1), 'r',phi_circl_rad, abs(j2), 'b');
            
            
            % ������� ����������� ��������� �� ��-� ������
            obj.j1 = j1;
            obj.j2 = j2;
            obj.Z = Z;
            obj.E = E;
            
        end
        
        function [H_calc] = H(obj, n, x)
            H_calc = besselh(n,1,x); 
        end

        function [J_calc] = J(obj, n,x)
            J_calc = besselJ(n,x); 
        end

        function [dJ_calc] = dJ(obj, n,x)
            dJ_calc = n/x*besselj(n,x) - besselj(n+1,x);
        end

        function [dH_calc] = dH(obj, n,x)
            dH_calc = n/x*besselh(n,1,x) - besselh(n+1,1,x);
        end
% %%        
%         % �������� ���� � ����
%         function obj = write_I_to_file(obj)
%             
%             % ������� 
%             I_real = obj.I_real;
%             N = obj.N;
%             n_for_graf = obj.n_for_graf;
%             
%             name_file = strcat('I_HFIE_MIE_Epol_sqrt_', num2str(N),'.dat');
%             
%             f01 = fopen(name_file,'w');
%             for n=1:N
%                 fprintf(f01,' %10.5f %10.5f\n', n_for_graf(n)*dx, abs(I_real(n)));
%             end         
%             fclose(f01);         
%         end
%%        
        % ��������� ���
        function obj= calculate_DA(obj)
            
            % ������� ����������� ��������� �� ��-� ������ 
                        % ������� ��������� �� ����
            x_midle = obj.circl_main.x_midle;
            z_midle = obj.circl_main.z_midle;
            dx = obj.dx;
            k2 = obj.k2;
            eta2 = obj.eta2;
            gamma = obj.gamma;
            N = obj.N;
            j2 = obj.j2;
            
            
            for p = 1 : 721
                Sum_E = 0;
                phi = (p-1)/2; % ���� �� ������� ��� ���������� 721 � 360 ��������
                phi_for_graf_DA(p) = (p-1)/2; % ���� ��� �������

                cx = - 1i*k2*cosd(phi);
                cz = - 1i*k2*sind(phi);


                    for n = 1 : N
                        xn = x_midle(n);
                        zn = z_midle(n);

                        Sum_E = Sum_E + eta2*k2/4*dx*j2(n)*exp(cx*xn+cz*zn);      
                    end

            RCS(p) = 10*log10((4/k2)*Sum_E*conj(Sum_E));

            end
            % ������ ���� � ������� ���� �� ���� �������� (�� �������������)
            figure;
            plot(phi_for_graf_DA, RCS,'r-');
            title('������ ������������� ���� � ������� ����'); 
            xlim([0 360])
            xlabel('����, ����'); 
            ylabel('���');
            
            % �������
            obj.DA = RCS; 
            obj.phi_for_graf_DA = phi_for_graf_DA;
%             obj.Sum_E_mass = Sum_E_mass;
        end
%  %%       
%         % ��������� ���� � ����
%         function obj= calculate_DE(obj)
%         end
%  %%       
%         % ��������� �������
%         function obj =  write_DA_to_file(obj)
%             
%             % �������
%             phi_for_graf_DA = obj.phi_for_graf_DA;
%             DA = obj.DA;
%             
%             f02 = fopen('DA_MFIE_MIE_Epol_sqrt.dat','w');
%             for p = 1 : 721
%                 fprintf(f02, ' %10.5f %10.5f\n', phi_for_graf_DA(p), DA(p));
%             end
%             fclose(f02);
%         end
%                 %%
%         % �������� ��� -  GETer
%         function I = get_I(obj)
%             I = obj.I_real;
%         end
%         % �������� ��� -  GETer
%         function DA = get_DA(obj)
%             DA = obj.DA;
%         end
    end
end