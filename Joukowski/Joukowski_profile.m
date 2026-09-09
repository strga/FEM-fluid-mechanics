%clear all
%clc

%Imaginary number
i = 1i;

%Parametres
a = 6;          %Real part
h = 0.5;        %Imaginary part 
delta = 0.6;    %Distance between centres of cyrcles M & K

U = 50; %Speed 
%gam = 0.5; %Cirkulace rychlosti po profilu
degUhel = 5; %Angle of attack - degrees

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Radius of circle Kr
R1 = sqrt(a.^2 + h.^2);

%Change of triangles (from similar triangles)
dh = delta * h / R1;
da = delta * a / R1;

Sx = - da; % centre of circle Mx ( minus profile orientation )
Sy = h + dh; % centre of circle My
R = R1 + delta; % Radius of M

%Joukowski profil
%tt = 2 * pi * ( 0 : 0.001 : 1 ); %0.001
tt = 2 * pi * ( 0 : 0.0104 : 1 );

z = Sx + R .* cos(tt) + i * ( Sy + R .* sin(tt) );

fceZ = (z + a.^2 ./ z) / 2;

figure(1)
plot( real(fceZ(1,:)) , imag(fceZ(1,:)) )
axis equal
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uhel = pi * degUhel / 180; %Angle of attack - radians THETA
Radius = (R:0.1:30)';

%Speed
Uinfad = U .* (cos(uhel) - sin(uhel) * i);
Uinf = U .* (cos(uhel) + sin(uhel) * i);

%Circlulation of velocity around the profile
gam = 2 * pi * R * U * ( uhel + atan(h/a) );

%Area in polar coordinates
[rr,phi] = meshgrid(Radius, tt );

%Polar complex number
xx = Sx + rr .* cos(phi);
yy = Sy + rr .* sin(phi);
z1 = (xx + i .* yy);

%Zukovskeho funkce -- body pro contour
zz = ( z1 + a.^2 ./ z1 ) / 2;

%zz_der = ( 1 - a.^2 ./ zz.^2 ) / 2;
%zz_der = ( 1 + ( z1 ./ ( sqrt( z1.^2 - a.^2 ) ) ) );
zz_der = 1/2 * ( 1 - ( a.^2 ./ z.^2 ) );

%zz2 = K(z) - S     K(z) = inverse mapping to Zuk. mapping
zz2 = z1 - Sx - i * Sy;

%Complex potential - Zukovskij
Fzz = ( 1/2 .* ( Uinfad * zz2 + Uinf * R.^2 ./ zz2 - gam / (pi * i) * log(zz2) ) );

%Complex speed - Zukovskij
Wzz = zz_der' .* ( 1/2 .* ( Uinfad - ( ( Uinf * R.^2 ) ./ ( zz2.^2 ) ) ) - gam / ( 2 * pi * i ) * ( 1 ./ zz2 ) );

tt = 0:0.0104:1; %Parametr 0.0099
aa = min(min(imag(Fzz))); %Minimal value of Im(Fzz) in the entire matrix
bb = max(max(imag(Fzz))); %Maximal value of Im(Fzz) in the entire matrix

%Plot streamlines
figure(2)
plot( real(fceZ(1,:)) , imag(fceZ(1,:)), 'color', 'k' );
hold on
contour(real(zz), imag(zz), imag(Fzz), aa + (bb - aa) * tt);
hold off
xlim([-8,7]); ylim([-7,7]);
%A = [real(fceZ(1,:)') , -imag(fceZ(1,:)'), zeros(1001,1)];

figure(3)
plot( real(fceZ(1,:)) , imag(fceZ(1,:)), 'color', 'k' );
hold on
contour(real(zz), imag(zz), sqrt( real(Wzz).^2 + imag(Wzz).^2 ), 100);
hold off
%xlim([-8,7]); ylim([-7,7]);
colorbar
