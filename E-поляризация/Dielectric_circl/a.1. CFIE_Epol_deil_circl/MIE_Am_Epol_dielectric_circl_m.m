%%  ���

close all;
clear variables;
%% 
% ������� ���������

%  1 - ������� ������ 
k0 = 2*pi;
% ��������� �����-1
e1 = 2.56;
mu1 = 1;
k1 = k0*sqrt(mu1*e1);
eta1 = 120*pi*sqrt(mu1/e1);

% ��������� �����-2
k2 = k0;
eta2 = 120*pi;

% ��������� ���������
a = 5/(2*pi);  % ������ �����
N = 200; % ���������� ����� �� �����
phi_i_grad = 0; % ���� ������� [����]
phi_i_rad = 2*pi/360*phi_i_grad; % ���� ������� [���]
gamma = 1.781072417990;

% ========================================================================
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
    phi_circl_rad(i) = (i-1)*dphi_rad + dphi_rad/2;
    phi_circl_grad(i) = (i-1)*dphi_grad + dphi_grad/2;
    a_for_graph(i) = a;
end
%% 
% ��������� �����

dphi_grad = 360/N;
x = zeros(1,N);
z = zeros(1,N);

for i = 1 : N
    z(i) = a*sind(dphi_grad*(i-1));
    x(i) = a*cosd(dphi_grad*(i-1));
end

plot(x,z)
%% 
% ��������� t ��� ������� �������

for i = 1 : N
    if i == N
        tx(i) = -((x(1)-x(i)))/dx;
        tz(i) = -((z(1)-z(i)))/dx;
    else
        tx(i) = -(x(i+1)-x(i))/dx;
        tz(i) = -(z(i+1)-z(i))/dx;
    end
end
%% 
% ��������� n ��� ������� �������

%  ��������� n ��� ������� �������
%  ���������� ������ ������ � ������� ������� ��������
    alpha = -90;
    matrix = [cosd(alpha),-sind(alpha); sind(alpha),cosd(alpha)];

for i = 1 : N
     matrix_temp = [tx(i), tz(i)]*matrix;
     nx(i) = matrix_temp(1);
     nz(i) = matrix_temp(2);
end
%% 
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

plot(x,z,x_midle, z_midle, 'rx');

% ������������� tx � tz
plot(tx, tz, 'rx');
compass(tx(1),tz(1));

% ========================================================================
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
            % c����������� E1y
            
            fm_y1 = +1/2;
            
            % c����������� E2y
            fm_y2 = -1/2;

            % c����������� H1u
            hm_u12 = 1i/4*dx*(1 + 1i*2/pi*(log(gamma*k1*dx/4)-1));
 
            % c����������� H2u
            hm_u22 = 1i/4*dx*(1 + 1i*2/pi*(log(gamma*k2*dx/4)-1));
            
        else
            % c����������� E1y
            fm_y1 = - 1i/4*k1*dx*(wmn/rmn*H(1,k1*rmn));
            
            % c����������� E2y
            fm_y2 = - 1i/4*k2*dx*(wmn/rmn*H(1,k2*rmn));

            % c����������� H1u
            hm_u12 = 1i/4*dx*H(0,k1*rmn);
            
            % c����������� H2u
            hm_u22 = 1i/4*dx*H(0,k2*rmn);
            
        end
        
      % ��������� ����
        % c����������� H1u
        hm_u11 = - 1i/4*k1*(umn_plus/r_plus*H(1,k1*r_plus) - umn_minus/r_minus*H(1,k1*r_minus));
        
        % c����������� H2u
        hm_u21 = - 1i/4*k2*(umn_plus/r_plus*H(1,k2*r_plus) - umn_minus/r_minus*H(1,k2*r_minus));
        
        % c����������� H1w
        hm_w1 = - 1i/4*k1*(wmn/r_plus*H(1,k1*r_plus) - wmn/r_minus*H(1,k1*r_minus));
        
        % c����������� H2w
        hm_w2 = - 1i/4*k2*(wmn/r_plus*H(1,k2*r_plus) - wmn/r_minus*H(1,k2*r_minus));
               
%%
        %% ������� ��� ������
%         fm_y1 = fm_y1;
%         fm_y2 = fm_y2;
%         hm_u11 = hm_u11;
%         hm_u21 = hm_u21;
%         hm_u12 = hm_u12;
%         hm_u22 = hm_u22;
%         hm_w1 = hm_w1;
%         hm_w2 = hm_w2;
        
        % c����������� E1y
        E1y = - fm_y1;
        
        % c����������� E2y
        E2y = - fm_y2;
       
        % c����������� H1u
        H1u = 1i/(eta1*k1)*(hm_u11 + (k1^2)*hm_u12);
        
        % c����������� H2u
        H2u = 1i/(eta2*k2)*(hm_u21 + (k2^2)*hm_u22);
        
        % c����������� H1w
        H1w = 1i/(eta1*k1)*hm_w1;
        
        % c����������� H2w
        H2w = 1i/(eta2*k2)*hm_w2;
  
%%
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
%%
% �������� ��� ������� 
plot(phi_circl_grad,  abs(j1), 'r',phi_circl_grad, abs(j2), 'b');

% �������� ������ �����
polarplot(phi_circl_rad, abs(j1), 'r',phi_circl_rad, abs(j2), 'b');
%% ���

for p = 1 : 721
    Sum_E = 0;
    phi = (p-1)/2; % ���� �� ������� ��� ���������� 721 � 360 ��������
    phi_for_graf_DA(p) = (p-1)/2; % ���� ��� �������
    
    cx = - 1i*k2*cosd(phi);
    cz = - 1i*k2*sind(phi);
    

        for n = 1 : N
            txn = tx(n);
            tzn = tz(n);
            xn = x_midle(n);
            zn = z_midle(n);
            
%             Sum_E = Sum_E + eta2*k2/4*dx*j2(n)*exp(cx*xn+cz*zn);
            Sum_E = Sum_E + k2*dx/4*j2(n)*(txn*sind(phi) - tzn*cosd(phi))*exp(cx*xn+cz*zn);

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

% ������������� ���� �� ����������� 
E_surf = zeros(N, 1); 
for m  = 1 : N
    xm = x_midle(m);
    zm = z_midle(m);
    txm = tx(m);
    tzm = tz(m);
    nxm = nx(m);
    nzm = nz(m);
    n_for_graf(m) = m; % ��� �������
    
    
    E_surf_n = 0;
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
           fm_y2 = -1/2;

        else
            % ������������� ����
            fm_y2 = - 1i/4*k2*dx*(wmn/rmn*H(1,k2*rmn));  
        end
        E2y = - fm_y2;
        E_surf_n = E_surf_n + j2(n) * E2y;
    end
    
    E_surf(m) =  E_surf_n; 
 end
 
figure;
plot(n_for_graf, abs(E_surf), 'r');
title('������������� ���� �� �����������'); 

% name_file = strcat('I_MIE_circl_diel_Epol_j1_', num2str(N),'.dat');
% f01 = fopen(name_file,'w');
% 
% % �� �������� Y ������ �������� � ������ �������� ������� �������
% 
% for n=1:N
%     fprintf(f01,' %10.5f %10.5f\n', phi_circl_grad(n), abs(j1(n)));
% end
% 
% fclose(f01);
%%
% name_file = strcat('I_MIE_circl_diel_Epol_j2_', num2str(N),'.dat');
% f02 = fopen(name_file,'w');
% 
% % �� �������� Y ������ �������� � ������ �������� ������� �������
% 
% for n=1:N
%     fprintf(f02,' %10.5f %10.5f\n', phi_circl_grad(n), abs(j2(n)));
% end
% 
% fclose(f02);
%%
% name_file = strcat('DA_MIE_circl_diel_Epol_j2_', num2str(N),'.dat');
% f03 = fopen(name_file,'w');
% 
% % �� �������� Y ������ �������� � ������ �������� ������� �������
% 
% for n= 1 : 721
%     fprintf(f03,' %10.5f %10.5f\n', phi_for_graf_DA(n), RCS(n));
% end
% 
% fclose(f03);
%%
function [H_calc] = H(n,x)
    H_calc = besselh(n,1,x); 
end

function [J_calc] = J(n,x)
    J_calc = besselJ(n,x); 
end
    
function [dJ_calc] = dJ(n,x)
    dJ_calc = n/x*besselj(n,x) - besselj(n+1,x);
end

function [dH_calc] = dH(n,x)
    dH_calc = n/x*besselh(n,1,x) - besselh(n+1,1,x);
end