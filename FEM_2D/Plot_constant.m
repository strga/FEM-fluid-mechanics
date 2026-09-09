function [Constants,Psi_point] = Plot_constant(M,bvp)

x = M.POS(:,1); y = M.POS(:, 2);
M.tri = M.TRIANGLES(1:M.nbTriangles,1:3);

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Constants = linspace(-10,10,20); %Values of constant for Dirichlet BC
bvp.MarkDir1 = 10;
Phys_point = find(M.POINTS(:,2) == 3333);
Point = M.POINTS(Phys_point,1); %Find point with physical value 3333
n = length(Constants);
Psi_point = zeros(n,1);

for i = 1:n

    j = Constants(i);
    bvp.Dir1 = @(x,y) j;
    
    [A,b] = DiscretizePotentialFlow(M,bvp);
    psi = A \ b;
    
    Psi_point(i) = psi(Point); %Hodnoty proudove funkce v bode v uzlu Point
end


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


end