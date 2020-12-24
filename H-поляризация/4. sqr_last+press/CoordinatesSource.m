function [ s ] = CoordinatesSource()

global as dx N

% Сооздадим структура для хранения координат источников X и Z
%s - source
s = struct('x', {[]}, 'z', {[]});

% генератора координат для первого участка
x1 = zeros (1,N);
z1 = zeros (1,N);
x2 = zeros (1,N);
z2 = zeros (1,N);
x3 = zeros (1,N);
z3 = zeros (1,N);
x4 = zeros (1,N);
z4 = zeros (1,N);

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

s(1).x = x1;
s(1).z = z1;
s(2).x = x2;
s(2).z = z2;
s(3).x = fliplr(x3);
s(3).z = z3;
s(4).x = x4;
s(4).z = fliplr(z4);

%временная переменная доя гарфика
temp_x = [s(1).x s(2).x s(3).x s(4).x];
temp_z = [s(1).z s(2).z s(3).z s(4).z];


% график точки невязок
figure;
% subplot(2,2,4);
plot(temp_x, temp_z,'-bo');
title('ГРАФИК_ФУНКЦИИ:CoordinatesSource - Точки расположение ВИ'); 

% ASSS = 15;

% clearvars ('x1', 'x2'); 

end


