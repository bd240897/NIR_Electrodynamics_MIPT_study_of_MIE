classdef D_class_MIE_MFIE_sqr_E_pol < P_class_MIE_sqr
%% """����� ��� MIE EFIE E-pol sqr"""
% 02-03-21

% ������ �������������:
    % as = 1;
    % phi_grad = 60;
    % N= 100;
    %     
    % calc_1 = D_class_MIE_MFIE_sqr_E_pol(as, phi_grad);
    % calc_1.calculate_main_sqr(N)
    % calc_1.calculate_I(N).write_I_to_file();
    % calc_1.calculate_DA().write_DA_to_file();
    % I = calc_1.get_I();
    % DA = calc_1.get_DA();
    
%% *** �������� ������ ***
%%  ��� ����������� ��������    
    %%  ��� ����������� ��������
    properties 

    end   
    
%% *** ������***
    methods  
%% ����������� - ������������
        % �����������
        function obj = D_class_MIE_MFIE_sqr_E_pol(as, phi_grad)
            obj@P_class_MIE_sqr(as, phi_grad);
            obj.type_MIE = string("MFIE");
            obj.type_pol = string("Epol");
            obj.name_of_class = strcat("class\_MFIE\_Epol\_sqr");
        end
%% ������ �����       
        % ������ ����� ��� ������� ����� ���������� � ����� �� ������
        function obj = calculate_I(obj, N)

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
            type_MIE = obj.type_MIE;
            type_pol = obj.type_pol;
            
            % �������� ������
            obj.print_sqr()
                               
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
                    umn = xmO2;
                    wmn = zmO2;                    

            %       ��������� ����� m � n (� ���)
                    xmn = xm - xn;
                    zmn = zm - zn;
                    r = sqrt(xmn^2 + zmn^2);

            %       ������ ��������� ��������� (� ���)
                    Hus = 0;
                    Hws = 0;
                    
                    if n == m %r < 0.001*dx 

                        % u - c����������� 
                        Hus = 0.5;

                        % w - c�����������
                        Hws = 0;

                    elseif abs(n - m) < 1.5 %r <= 2*dx
                        % u - c�����������                     
                        umn_plus_dx2 = umn + dx/2;
                        umn_minus_dx2 = umn - dx/2;
                        Hus = (1i/8*wmn*k^2*dx + 1/(2*pi)*(atan(umn_plus_dx2/wmn) - atan(umn_minus_dx2/wmn)));

                        % w - c�����������
                        r_plus_dx2 = sqrt((umn + dx/2)^2 + wmn^2);
                        r_minus_dx2 = sqrt((umn - dx/2)^2 + wmn^2);
                        Hws = 1i/4*(besselh(0,1,k*r_plus_dx2) - besselh(0,1,k*r_minus_dx2));

                    else
                        % u - c�����������
                        rmn = sqrt(umn^2 + wmn^2);    
                        Hus = 1i/4*k*dx*(wmn/rmn*besselh(1,1,k*rmn));  

                        % w - c�����������
                        r_plus_dx2 = sqrt((umn + dx/2)^2 + wmn^2);
                        r_minus_dx2 = sqrt((umn - dx/2)^2 + wmn^2);
                        Hws = 1i/4*(besselh(0,1,k*r_plus_dx2) - besselh(0,1,k*r_minus_dx2));
                    end

            %       ������ ��������
                    if m == n
                        kronker = 1;
                    else 
                        kronker = 0;
                    end

                    % ������ �������� ��������� 
                    Zmn(m,n) = kronker + Hus*(nxm*tzn - nzm*txn) + Hws*(nxm*nzn - nzm*nxn);

                end
            end

            % ������ ��������� ����
            for m= 1 : 4*N
                nxm = nx(m);
                nzm = nz(m);

                cx = -1i*k*cosd(phi_grad);
                cz = -1i*k*sind(phi_grad);

                % ��������� ����������
                xm = X(m);
                zm = Z(m);

            %   �������� �������� ���� - ��� ������������ �����
            %     etta = 120*pi;
                Hxi = 1*sind(phi_grad)*exp(cx*xm+cz*zm);
                Hzi = -1*cosd(phi_grad)*exp(cx*xm+cz*zm);

                Hi(m,1) = -nxm*Hzi + nzm*Hxi;% �������� ���� X
            end

            % ��������� ���� � �������� ������
            % �������� ����� �� 
            I=Zmn\Hi;
            
            % ������ �� ��� ������� 
            n_for_graf = zeros(1, N);
            for i = 1 : 4*N
                 n_for_graf(i) = i;
            end
            
            % ��� X � ������������� ��������
            figure;
            plot(n_for_graf*dx, abs(I), 'r');
            

            
            
            title("������ ������������� ���� (��� ����� 2)  - (4)"); 
            xlabel("��������"); 
            ylabel("��������� ����");
            legend(obj.name_of_class);
            
            % ������� ����������� ��������� �� ��-� ������
            obj.I_real = I; 
            obj.n_for_graf = n_for_graf;
            obj.Zmn = Zmn;
            obj.Hi = Hi;
            
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
            type_MIE = obj.type_MIE;
            type_pol = obj.type_pol;
            
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
            title("������ ������������� ���� � ������� ����"); 
            xlim([0 360])
            xlabel("����, ����"); 
            ylabel("���");
            legend(obj.name_of_class);
            
            % �������
            obj.DA = RCS; 
            obj.phi_for_graf_DA = phi_for_graf_DA;
            obj.Sum_E_mass = Sum_E_mass;
           
        end
    end    
end