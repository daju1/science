function [x,y,z] = generate_random_sphere_points(npoints, r)

fi = 2*pi*rand([1,npoints]);
theta = acos(1-2*rand([1,npoints]));

x = r * sin(theta) .* cos(fi);
y = r * sin(theta) .* sin(fi);
z = r * cos(theta);

end