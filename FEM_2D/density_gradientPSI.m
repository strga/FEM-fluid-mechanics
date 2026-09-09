function [U1, rho_grad] = density_gradientPSI(max_speed, rho)

U1 = 0:0.5:max_speed; %Speed/gradient of PSI
n = length(U1);
%rho = 1.225; %initial value of density

for k = 1:n
    U = [U1(k), 0];
    rho = zavislost(U, rho); % Calculates value of density by itereation for given speed U [ grad(PSI) ]
    rho_grad(k) = rho;
end

end
%plot(U1, rhog);