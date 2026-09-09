clear all;
clc;

i = 1i; %imaginary number

R = 0.1; %Radius of cylinder
gam = 0.05; %Circulation of velocity around profile
degUhel = 8; %Angle of attack - degrees

uhel = pi * degUhel / 180; %Angle of attack - radians

%Speed
Uinfad = 10 * (cos(uhel) - sin(uhel) * i);
Uinf = 10 * (cos(uhel) + sin(uhel) * i);

%Region Omega - polar coordinates, radius of contour area r = < 0.1 ; 1.1 >
[rr,phi] = meshgrid(R + 10 * R * (0:0.002:1), 0:0.001:(2 * pi) );

%Polar complex number
xx = rr .* cos(phi);
yy = rr .* sin(phi);
zz = xx + i .* yy;

%Complex potential
Fzz = Uinfad * zz + Uinf * R^2 ./ zz - gam * log(zz) ./ (2 * pi * i);

tt = 0:0.02:1; %Parametr
aaa = min(min(imag(Fzz))); %Minimal value of Im(Fzz) in the entire matrix
bbb = max(max(imag(Fzz))); %Maximal value of Im(Fzz) in the entire matrix

%Plot -> xx, yy, Im(Fzz) = proudova fce, levels from min. to max. value
contour(xx, yy, imag(Fzz), aaa + (bbb - aaa) * tt);
