function rho = zavislost(U,rhoinit)

rho = rhoinit;

err = 1;

while (err > 1e-12)
    irho = invrho(U,rho); %Calculation of (1 / rho) from Bernoulli density equation, speed U == grad( PSI )
    rho1 = 1. / irho;
    err = abs(rho - rho1);
    rho = rho1;
end