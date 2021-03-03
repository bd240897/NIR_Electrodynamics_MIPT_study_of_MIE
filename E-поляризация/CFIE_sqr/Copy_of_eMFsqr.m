% Scattering of a plane wave by a perfectly conducting cylinder of 
% SQUARE cross section.
% E-polarization.
% Method of Electric Field Intergral Equation (EFIE);
% Method of Magnetc Field Intergral Equation (MFIE);
% Method of Combined Field Intergral Equation (CFIE).

%f01 = fopen('geo.dat','w');
%f03 = fopen('hdE1n.dat','w');
%f05 = fopen('hRCS1nb.dat','w');
f07 = fopen('aI.dat','w');
tic;

k  = 2*pi;
gk2 = 1.781072*k/2;
et0 = 120*pi;
aph = 0.5;
ap1 = 1-aph;
a  = 2; %/pi; %0.3;     % width of square in wavelengths
ka = k*a;
%for ni = 1:361
phi = 0;     % angle of incidence from the x axis, degrees
phr = pi*(phi/180);
sph = sin(phr);
cph = cos(phr);
ikc =-1i*k*cph;
iks =-1i*k*sph;

N1 = 100; % number of segments on one side
N2 = 2*N1;
N3 = 3*N1;
N  = 4*N1;

x = zeros(N,1);
z = x;
s = x;
ds = x;
tx = x;
tz = x;
nx = x;
nz = x;

d  = a/N1;
hd = d/2;
for n = 1:N1
    nd = (n-0.5)*d;
    n2 = n+N1;
    n3 = n+N2;
    n4 = n+N3;
    x(n) = 0;
    z(n) = a-nd;
    tx(n) = 0;
    tz(n) = 1;
    nx(n) =-1;
    nz(n) = 0;
    ds(n) = d;
    s(n)  = nd;
    x(n2) = nd;
    z(n2) = 0;
    tx(n2) =-1;
    tz(n2) = 0;
    nx(n2) = 0;
    nz(n2) =-1;
    ds(n2) = d;
    s(n2)  = nd+a;
    x(n3) = a;
    z(n3) = nd;
    tx(n3) = 0;
    tz(n3) =-1;
    nx(n3) = 1;
    nz(n3) = 0;
    ds(n3) = d;
    s(n3)  = nd+2*a;
    x(n4) = a-nd;
    z(n4) = a;
    tx(n4) = 1;
    tz(n4) = 0;
    nx(n4) = 0;
    nz(n4) = 1;
    ds(n4) = d;
    s(n4)  = nd+3*a;
    %fprintf(f01,' %10.5f %10.5f\n',x,z);
end

% System of algebraic equations
i4 = 1i/4;
ZM = zeros(N,N);
H  = zeros(N,1);
ZE = zeros(N,N);
E  = zeros(N,1);
for m  = 1:N
    xm = x(m);
    zm = z(m);
    ZM(m,m) = 0.5;
    ZE(m,m) = (k*ds(m)/4)*(1+(2i/pi)*(log(gk2*ds(m)/2)-1));
    for n  = 1:N
        if abs(m-n) > 0.5
            bu = nz(m)*tx(n)-nx(m)*tz(n);
            bw = nz(m)*nx(n)-nx(m)*nz(n);
            dx = xm-x(n);
            dz = zm-z(n);
            u  = dx*tx(n)+dz*tz(n);
            w  = dx*nx(n)+dz*nz(n);
            w2 = w^2;
            u1 = u+ds(n)/2;
            u2 = u-ds(n)/2;
            r1 = sqrt(u1^2+w2);
            r2 = sqrt(u2^2+w2);
            fw = i4*(besselh(0,1,k*r2)-besselh(0,1,k*r1));
            if abs(n-m) > 1.5
                r  = sqrt(u^2+w2);
                kr = k*r;
                fu = k*ds(n)*(w/r)*besselh(1,1,kr);
                ZE(m,n) = (k*ds(n)/4)*besselh(0,1,kr);
            else
                fu = (2i/pi)*(atan2(u2,w)-atan2(u1,w));
                fu = k*w*k*ds(n)/2+fu;
                fe = u1*log(gk2*r1)-u2*log(gk2*r2)-ds(n);
                fe = fe+w*(atan2(u1,w)-atan2(u2,w));
                ZE(m,n) = (k/4)*(ds(n)+(2i/pi)*fe);
            end
            fu =-i4*fu;
            ZM(m,n) = bu*fu+bw*fw;
        end
    end
    g = exp(ikc*xm+iks*zm);
    H(m,1) = (nz(m)*sph+nx(m)*cph)*g;
    E(m,1) = g;
end
toc;
tic;

IM = ZM\H;
IE = ZE\E;
ZC = aph*ZE+ap1*ZM;
C  = aph*E +ap1*H;
IC = ZC\C;
toc;

for n  = 1:N
    aIe = abs(IE(n,1));
    aIm = abs(IM(n,1));
    aIc = abs(IC(n,1));
    xx  = s(n)/(4*a);
    fprintf(f07,' %10.5f %10.5f %10.5f %10.5f\n',xx,aIm,aIe,aIc);
end

%Residual

% Radar cross-section
for np = 1:1%721
    pd = phi;%(np-1)/2;
    ph = pi*(pd/180);
    cop = cos(ph);
    sip = sin(ph);
    f = 0;
    for n = 1:N
        f = f;
    end
    af = abs(f*sip);
    sigma = 2*af*af/pi;
    if sigma < 0.001
        sigma = 0.001;
    end
    sigma = 10*log10(sigma);
    %fprintf(f05,' %10.3f %12.6f %10.3f\n',pd,sigma,phi);
end

%end

%fclose(f01);
%fclose(f03);
%fclose(f05);
fclose(f07);
