clear;clc;

% скорость распространения потенциала
cd = 300000000;

% шаг времени
dt = 1e-15;
vt = 0:dt:1e-9;

% История координат заряда
x = 150000000*vt;
y = zeros(size(vt));
z = zeros(size(vt));

% История скоростей заряда
vx = zeros(size(vt));
vy = zeros(size(vt));
vz = zeros(size(vt));

FI = [];
XFI = [];
YFI = [];
ZFI = [];

% Рассчитываем картину распределения запаздывающего потенциала в момент времени 
i_t = size(vt, 2);
t = vt(i_t);

dr = cd * dt;

% В цикле по истории движения заряда
vi = flip(2:i_t-1);
for i = vi
    % координата заряда в момент "излучения" потенциала
    xe = x(i);
    ye = y(i);
    ze = z(i);
    
    % скорость заряда в момент "излучения" потенциала
    vx(i) = (x(i+1)-x(i-1))/(2*dt);
    vy(i) = (y(i+1)-y(i-1))/(2*dt);
    vz(i) = (z(i+1)-z(i-1))/(2*dt);
    
    % время прошедшее с момента "излучения" потенциала до момента
    % построения картины распределения
    T = t - vt(i);
    
    % расстояние на которое произошло распространение потенциала с момента
    % его "излучения" до момента построения картины распределения
    R = cd * T;
    
    
    % отсекаем не интересующую нас область пространства
    if R > 300
        return
    end
    
    npoints = 800;
    
    [xfi, yfi, zfi, fi] = generate_mende_potencial_points(npoints, R, xe, ye, ze, vx(i), vy(i), vz(i));
    
    XFI = [XFI, xfi];    
    YFI = [YFI, yfi];    
    ZFI = [ZFI, zfi];
    
    FI = [FI, fi];
    
end

vx(1) = (-3*x(1)+4*x(1+1)-x(1+2))/(2*dt);
vy(1) = (-3*y(1)+4*y(1+1)-y(1+2))/(2*dt);
vz(1) = (-3*z(1)+4*z(1+1)-z(1+2))/(2*dt);

vx(i_t) = (3*x(i_t)-4*x(i_t-1)+x(i_t-2))/(2*dt);
vy(i_t) = (3*y(i_t)-4*y(i_t-1)+y(i_t-2))/(2*dt);
vz(i_t) = (3*z(i_t)-4*z(i_t-1)+z(i_t-2))/(2*dt);


vtt = vt(i_t-1) : dt/5 : vt(i_t);
i_tt = size(vtt, 2);


x_tt = linterp(vt, x, vtt);
y_tt = linterp(vt, y, vtt);
z_tt = linterp(vt, z, vtt);

vx_tt = linterp(vt, vx, vtt);
vy_tt = linterp(vt, vy, vtt);
vz_tt = linterp(vt, vz, vtt);



% В цикле по истории движения заряда
vi_tt = flip(2:i_tt-1);
for i = vi_tt
    % координата заряда в момент "излучения" потенциала
    xe = x_tt(i);
    ye = y_tt(i);
    ze = z_tt(i);
    
    % время прошедшее с момента "излучения" потенциала до момента
    % построения картины распределения
    T = t - vtt(i);
    
    % расстояние на которое произошло распространение потенциала с момента
    % его "излучения" до момента построения картины распределения
    R = cd * T;
    
    
    % отсекаем не интересующую нас область пространства
    if R > 300
        return
    end
    
    npoints = 1500;
    
    [xfi, yfi, zfi, fi] = generate_mende_potencial_points(npoints, R, xe, ye, ze, vx_tt(i), vy_tt(i), vz_tt(i));
    
    XFI = [XFI, xfi];    
    YFI = [YFI, yfi];    
    ZFI = [ZFI, zfi];
    
    FI = [FI, fi];
    
end


close all

dd  = dr/5;
l1 = 0.35;
l2 = 0.15;
% dx = floor(min(XFI)):dd:ceil(max(XFI));
% dy = floor(min(YFI)):dd:ceil(max(YFI));
% dz = floor(min(ZFI)):dd:ceil(max(ZFI));

dx = x(end) -l1 :dd : x(end)+l2;
dy = y(end) -l1 :dd : y(end)+l1;
dz = z(end) -l1 :dd : z(end)+l1;

[xq,yq,zq] = meshgrid(dx, dy, dz);

fiq = griddata(XFI, YFI, ZFI, FI, xq,yq,zq);

figure(1)
xslice = []; yslice = 0; zslice = [0];
slice(xq,yq,zq, fiq, xslice,yslice,zslice);
colormap hsv

% figure (2)
hold on
plot3(x,y,z)

figure(3)
contourslice(xq,yq,zq, fiq, xslice,yslice,zslice);
figure(4)
patch(isosurface(xq, yq, zq, fiq))
% figure(5)
% scatter3(XFI, YFI, ZFI, FI)



