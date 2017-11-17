function [xfi, yfi, zfi, fi] = generate_mende_potencial_points(npoints, r, xe, ye, ze, vxe, vye, vze)  

% ���������� ������ � ������ "���������" ����������
% xe ye ze

% �������� ������ � ������ "���������" ����������
% vxe vye vze

% ���������� �� ������� ��������� ��������������� ���������� � �������
% ��� "���������" �� ������� ���������� ������� �������������
% r


% �������� �����
c = 300000000;

% ���������� ������ ����� ��������� �� ����� "���������" ���������� �� ���������� ��� ���������������
%     [X,Y,Z] = sphere(100);
%     
%     rx = R * X(:);
%     ry = R * Y(:);
%     rz = R * Z(:);
%     

[rx,ry,rz] = generate_random_sphere_points(npoints, r);

% ������� ������� ���������� ������������ ������� �������� �� ������ ������
% syms i j k vxe vye vze rx ry rz
% det ([i j k; vxe vye vze; rx ry rz])
% 
% ans =
%  
% i*rz*vye - i*ry*vze + j*rx*vze - j*rz*vxe - k*rx*vye + k*ry*vxe

VxR_i2 = (rz*vye - ry*vze).^2;
VxR_j2 = (rx*vze - rz*vxe).^2;
VxR_k2 = (rx*vye + ry*vxe).^2;

fi = (1/r) * cosh(sqrt(VxR_i2 + VxR_j2 + VxR_k2) / ( r * c));

xfi = xe + rx;
yfi = ye + ry;
zfi = ze + rz;