clc;
clear all;
format long

n1 = 10000;     %kroky pro 'analyticke' reseni
n2 = 5;         %kroky pro numericke reseni MKP
K = [2 4 8 16 32 64 128 256 512];

f = @(x) x + 2 ;

%'Analyticke' reseni
[x,Y] = FEM_1D(n1,f);

%Numericke reseni
[a,b] = FEM_1D(n2,f); 

%'Analyticke' reseni - H = x-value, Int = y-value 
H = (0:(1/n1):1)';
Int = interp1(a,b,H);

xx = (1:1:n1+1)'; % x-value of 'analytic' solution

Q = Int(xx); %y-value of 'analytic' solution

Chyba = Y - Q; %y_numerical - y_'exact'
Max_chyba = max(Chyba);
  
chyba_h = [];
vel_h = [];

for n3 = K
    
    [r,s] = FEM_1D(n3,f); %vysledky [x,y]
    Int3 = interp1(r,s,H); %Int3 = Y-value, H = x-value
    GG = max( abs( Y - Int3(xx) ) ); % max error for given step

    chyba_h = [chyba_h, GG]; %souradnice y
    vel_h = [vel_h, 1. / (n3)]; %souradnice x

end

%Vypocet radu chyby p
for i = 2:length(K)
   
   p(1) = 0;
   j = chyba_h(i) ./ chyba_h(i-1);
   p(i) = -log2(j);
   
end

Krok_chyba = [vel_h; chyba_h; p]';

% Srovnani analytickeho reseni s MKP
figure(1)
plot(x,Y)
hold on
plot(a,b)
%title('Porovnání analytického řešení s MKP')
legend('Analytické řešení pro f = x + 2',"MKP s krokem h = 1/" + n2)

% Velikost chyby metody MKP proti analytickemu reseni
figure(2)
plot(H,Chyba,'-','MarkerSize',15)
%title("Chyba výpočtu MKP s krokem h = 1/" + n2)


% Velikost chyby v zavislosti na velikosti kroku h
figure(3)
loglog(Krok_chyba(:,1), Krok_chyba(:,2),'.','MarkerSize',15)
%figure(4)
%plot(Krok_chyba(:,1), Krok_chyba(:,2),'.','MarkerSize',15)
%title("Chyba výpočtu MKP podle velikosti kroku h v log-log souřadnicích s krokem h = 2^-n pro n = 1..." + K)
xlabel("Krok h / log")
ylabel("Velikost chyby E / log")
hold on

%Calculate slope M of log-log graph

M = log( chyba_h(3) / chyba_h(4) ) / log( vel_h(3) / vel_h(4) );
MM = log( chyba_h(6) / chyba_h(7) ) / log( vel_h(6) / vel_h(7) );

P = chyba_h(4) / chyba_h(3);
PP = chyba_h(3) / chyba_h(2);

% Porovnání chyby MKP
function [x,Y] = FEM_1D(n1, f)

% Rozdělení intervalu

n = n1 - 1;             % rozdělení na množství hodnot, které chci zjistit ... dim V = n
h = 1. / (n1);          % velikost kroku ... množství kroků, potřebné dostat se z v(0) do v(1), když hledám n hodnot mezi v(0) a v(1) (= dim V)
x = (0:h:1)';           % diskretizace x po velikosti h od 0 do 1
y = f(x);               % výsledné hodnoty z funkce f = @(x) 

% Matice tuhosti

U = diag(2 * ones(n,1));
V = diag(-ones(n-1,1),1);
W = diag(-ones(n-1,1),-1);
A = U + V + W;                  % Řídká matice tuhosti 

b = h * h * f(x(2:end - 1));    % Výpočet elementů matice od druhé do předposlední hodnoty (1. a poslední hodnota dána okrajovým podmínkami).
                                % Matice má tvar pro h = 1/(1+n) 1/h Au = b a zároveň b = h*f(x)'
u = A \ b;

Y = [0;u;0];                    % Vektor sloupcový s prvním elementem 0, pak elementy u a poslední 0

end


%{

clear all;
n1 = 10000;
n2 = 5;
f = @(x) x + 2;
[x,Y] = FEM_1D(n1,f);
[a,b] = FEM_1D(n2,f);

H = (0:(1/n1):1)';
Int = interp1(a,b,H);

xx = (1:1:n1+1)';

Q = Int(xx);  %Q je redundantní kvůli povaze Int
Chyba = Y - Q;
Max_chyba = max(Chyba);
  
chyba_h = [];
vel_h = [];

for n3 = 1:100
    
    [r,s] = FEM_1D(n3,f);
    Int3 = interp1(r,s,H);
    GG = max(Y - Int3(xx));

    chyba_h = [chyba_h, GG];
    vel_h = [vel_h, 1. / (n3)];

end

Krok_chyba = [vel_h; chyba_h]';

figure(1)

plot(x,Y)
hold on
plot(a,b)
title('Porovnání analytického řešení s MKP')
legend('Analytické řešení pro f = x + 2',"MKP s krokem h = 1/" + n2)

figure(2)

plot(H,Chyba,'-','MarkerSize',15)
title("Chyba výpočtu MKP s krokem h = 1/" + n2)

figure(3)

loglog(Krok_chyba(:,2), Krok_chyba(:,2))
title("Chyba výpočtu MKP podle velikosti kroku h v log-log souřadnicích")
xlabel("Krok h / log")
ylabel("Velikost chyby E / log")

%}