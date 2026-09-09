clear all;
clc;

i = 1i; %imaginary number

R = 0.1; %Radius of cylinder
gam = 0.05; %Cirkulacion of velicity around profile
degUhel = 8; %Angle of attack - degrees

uhel = pi * degUhel / 180; %(Angle of attack) - radians

%Speed
Uinfad = 10 * (cos(uhel) - sin(uhel) * i);
Uinf = 10 * (cos(uhel) + sin(uhel) * i);

%Region Omega
[xx,yy] = meshgrid( -3:0.01:3, -3:0.01:3 );

%Complex number
zz = xx + i * yy;

%Complex potential
Fzz = Uinfad * zz + Uinf * R^2 ./ zz - gam * log(zz) ./ (2 * pi * i);

% Complex potential inside cylinder (circle) area must be 0 = solid object
Fzz(xx.^2 + yy.^2 < R^2) = 0;

tt = 0:0.005:1; %Parametr
aaa = min(min(imag(Fzz))); %Minimal value of Im(Fzz) in the entire matrix
bbb = max(max(imag(Fzz))); %Maximal value of Im(Fzz) in the entire matrix

%Plot -> xx, yy, Im(Fzz) = proudova fce, levels from min. to max. value 
contour(xx, yy, imag(Fzz), aaa + (bbb - aaa) * tt);
