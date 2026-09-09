 function ret = invrho(U,rho)

kPa = 1e3;
rho0 = 1.225;
p0 = 101.325 * kPa;
gamma = 1.4; % 2 atom idealni plyn - Poisson constant
a0sqrt = gamma * p0 / rho0;
velsqrt = U(1) * U(1) + U(2) * U(2);

gam1 = gamma - 1;
ret = ( (1 - 1. / 2. * gam1 / a0sqrt * velsqrt ./ rho ./ rho).^gam1 / rho0 ); %Vztah pro inverse density z Beroulliho vzorce
end