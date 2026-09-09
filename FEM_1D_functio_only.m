% Porovnání chyby MKP
function [x,Y] = FEM_1D(n1, f)

%f = @(x) x + 2;

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