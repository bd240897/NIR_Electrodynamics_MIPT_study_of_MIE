classdef D_class_MIE_EFIE_sqr_E_pol < P_class_MIE_sqr
%% """����� ��� MIE EFIE E-pol sqr"""
% 02-03-21

% ������ �������������:
    % as = 1;
    % phi_grad = 60;
    % N= 100;
    %     
    % calc_1 = D_class_MIE_EFIE_sqr_E_pol(as, phi_grad);
    % calc_1.calculate_main_sqr(N)
    % calc_1.calculate_I(N).write_I_to_file();
    % calc_1.calculate_DA().write_DA_to_file();
    % I = calc_1.get_I();
    % DA = calc_1.get_DA();
    
%% *** �������� ������ ***
%%  ��� ����������� ��������    
    %%  ��� ����������� ��������
    properties 
        % ����� ���� ��� - ��� ���� � ������������ ������
    end   
    
%% *** ������***
    methods  
%% ����������� - ������������
        % �����������
        function obj = D_class_MIE_EFIE_sqr_E_pol(as, phi_grad)
            obj@P_class_MIE_sqr(as, phi_grad);
            obj.type_MIE = "EFIE";
            obj.type_pol = "Epol";
            obj.name_of_class = strcat("class\_EFIE\_Epol\_sqr");
        end
%% ��������� D0        

        % ��������� D0      
        function D0 = ganerate_D0(obj)
            % ��������� D0
            % �������
            dx = obj.dx;
            k = obj.k;
            gamma = obj.gamma;
            
            % ����� ��������
            kd4 = k*dx/4;
            D0 = k*dx*(1+(1i*2/pi)*(log(gamma*kd4)-1));
        end     
        
%% ������ �����       
        % ������ ����� ��� ������� ����� ���������� � ����� �� ������
        function obj = calculate_I(obj, N)
            
            % ������ ���� ����
            % ��������� ����� ��� �������
            obj.calculate_main_sqr(N);
            % ������� ��������� �� ����
            % ������� ��������� �� ����
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
            
            % �������� ������
            obj.print_sqr()
            
            % c���������� DO ��� ������� 
            D0 = ganerate_D0(obj);
            
            % ��������� ������ 
            delta = dx/2;
            
            % ������ ���� ����
            count = 0;
            for m = 1 : 4*N
            %   ���������� � �������� ������� ��������� (� ���) 
            % - �.�. ����� ��������� � ������ �������� ���������
                xm = X(m);
                zm = Z(m);
                txm = tx(m);
                tzm = tz(m);
                nxm = nx(m);
                nzm = nz(m);

                for n = 1 : 4*N
                %   ���������� � �������� ������� ��������� (� ���)
                    xn = X(n);
                    zn = Z(n);
                    txn = tx(n);
                    tzn = tz(n);
                    nxn = nx(n);
                    nzn = nz(n);

                %   ����� ����� ������� �������� - ������ ������� (� ���)
                %   ����� ������� �������� �������
                    x_begin_O2 = xn;
                    z_begin_O2 = zn;

                %   ��������� ����� ������� ��������� �2 � ������ ���������� (� ���)
                    x_mO2 = xm - x_begin_O2;
                    z_mO2 =  zm - z_begin_O2;

            %       ������� ���� ������ �� ����� ��
            %       ���������� ���� ����� (� ���)
                    xmO2 = x_mO2*txn + z_mO2*tzn;
                    zmO2 = x_mO2*nxn + z_mO2*nzn;

            %       �.�. ����� �� ������ �������� ������� �� ��� �����
                    % ���������� � ��������� �����������
%                     umn = xmO2;
%                     wmn = zmO2;                    
                    xmnO2 = xmO2;
                    zmnO2 = zmO2;

            %       ��������� ����� m � n (� ���)
                    xmn = xm - xn;
                    zmn = zm - zn;
                    r = sqrt(xmn^2 + zmn^2);

            %       ������ ��������� ��������� (� ���)
                    if n == m %r < 0.001*dx 
                        Es(m,n) = D0;

                    elseif abs(n - m) < 1.5 %r <= 2*dx
                        xmn_plus_delta = (xmnO2) + delta;
                        xmn_minus_delta = (xmnO2) - delta;
                        r_plus_delta = sqrt(xmn_plus_delta^2 + zmnO2^2);
                        r_minus_delta = sqrt(xmn_minus_delta^2 + zmnO2^2);
            %             Es(m,n) = k*dx+k*2*1i/pi*(dx*log(gamma*k/2)-dx*zmn*(atan(xmn_plus_delta/zmn)- atan(xmn_minus_delta/zmn)) + 0.5*xmn_plus_delta*log(r_plus_delta^2)  -  0.5*xmn_minus_delta*log(r_minus_delta^2)- dx);
                        Es(m,n) = k*dx+k*2*1i/pi*(zmnO2*(atan(xmn_plus_delta/zmnO2) - atan(xmn_minus_delta/zmnO2)) + xmn_plus_delta*log(gamma*k/2*r_plus_delta) - xmn_minus_delta*log(gamma*k/2*r_minus_delta) - dx);

                    else
                        f2 = besselh(0,1,k*r);
                        Es(m,n) = k*dx*f2;
                    end
                end
            end

            % ������ ��������� ����
            for m=1: 4*N
                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);
                
                % ��������� ����������
                xm = X(m);
                zm = Z(m);

            %   �������� �������� ����
                Ei(m,1) = +4*exp(cx*xm+cz*zm);% �������� ���� X
            end

            % ��������� ���� � �������� ������
            % �������� ����� �� 
            I=Es\Ei;
            
            % ������ �� ��� ������� 
            n_for_graf = zeros(1, N);
            for i = 1 : 4*N
                 n_for_graf(i) = i;
            end
            
            % ��� X � ������������� ��������
            figure;
            plot(n_for_graf*dx, abs(I), 'r');
            title('������ ������������� ���� (��� ����� 2)  - (4)'); 
            xlabel('��������'); 
            ylabel('��������� ����');
            legend(obj.name_of_class);
            
            % ������� ����������� ��������� �� ��-� ������
            obj.I_real = I; 
            obj.n_for_graf = n_for_graf;
            obj.Zmn = Es;
            obj.Ei = Ei;
           
        end
        
%% ������ ���       
        % ��������� ���
        function obj = calculate_DA(obj)
            
            % ������� ����������� ��������� �� ��-� ������ 
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
                phi = (p-1)/2; % ���� �� ������� ��� ���������� 721 � 360 ��������
                phi_for_graf_DA(p) = (p-1)/2; % ���� ��� �������

                cx = - 1i*k*cosd(phi);
                cz = - 1i*k*sind(phi);


                    for n = 1 : 4*N
                        txn = tx(n);
                        tzn = tz(n);
                        xn = X(n);
                        zn = Z(n);
            
                        % ������ �� �� 120 �� *(120*pi)
                        Sum_E = Sum_E + k/4*dx*I_real(n)*exp(cx*xn+cz*zn);      

                    end
                    
            Sum_E_mass(p) = Sum_E;
            RCS(p) = 10*log10((2/pi)*Sum_E*conj(Sum_E));

            end

            % ������ ���� � ������� ���� �� ���� �������� (�� �������������)
            figure;
            plot(phi_for_graf_DA ,RCS,'r-');
            title('������ ������������� ���� � ������� ����'); 
            xlim([0 360])
            xlabel('����, ����'); 
            ylabel('���');
            legend(obj.name_of_class);
            
            % �������
            obj.DA = RCS; 
            obj.phi_for_graf_DA = phi_for_graf_DA;
            obj.Sum_E_mass = Sum_E_mass;
            
        end
    end    
end