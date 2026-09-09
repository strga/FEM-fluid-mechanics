clear all
clc

%M = load_gmsh('squareDIR.msh');
%M = load_gmsh('square_2_Dir_Neu.msh');
%M = load_gmsh('squareDIR_NEU_1.msh');
%M = load_gmsh('squareDIR_ROB_1.msh');
%M = load_gmsh('squareDIR_NEU_ROB.msh');
M = load_gmsh('semicircleDIR.msh');
%M = load_gmsh('semicircleNEU.msh');
%M = load_gmsh('semicircleROB.msh');

x = M.POS(:,1); y = M.POS(:, 2);
M.tri = M.TRIANGLES(1:M.nbTriangles,1:3);

% figure(1)
% trisurf(M.tri, x, y, 2 .* pi.^2 .* sin(pi*x) .* sin(pi*y));
% title('Function f(x,y)') 

%RHS of Poisson equation
f = @(x, y) 2 * pi^2 * sin(pi * x) .* sin(pi * y);

%Dirichlet boundary condition
Dir1 = @(x, y) sin(pi * x) .* sin(pi * y);

%Neumann boundary condition
psi = @(x, y, n1, n2) ( n1 * pi * cos(pi * x) .* sin(pi * y) + n2 * pi * sin(pi * x) .* cos(pi * y) );

%Robin boundary condition
alfa = 1; 
ur = @(x, y, n1, n2) 1 * ( (( n1 * pi * cos(pi * x) .* sin(pi * y) + n2 * pi * sin(pi * x) .* cos(pi * y) ) / alfa) + sin(pi * x) .* sin(pi * y) );

k = 1;

Number_triangles = M.nbTriangles;
Parametre_triange = M.TRIANGLES(k,:);    % Informace o elementu k
Size_triangle = size(M.TRIANGLES);       % Velikost matice
M.TRIANGLES(1:10,:);                     % Vytáhnout informace o 1:n prvcích
M.LINES(1:10,:);                         % Vytáhnout inforamce o 1:n křivkách

Element_info = getMeshElement(M,k);

format long
Plocha = getDomainVolume(M);
Integral_teziste = computeIntegral1(M, f);
Integral_vrcholy = computeIntegral2(M, f);
Integral_strany = computeIntegral3(M, f);
Boundary_points = BoundaryPts(M);

[A,b] = Matice_tuhosti(M,f,Dir1,psi,alfa,ur);

%Numerical solution
u = A \ b;

figure(2)
trimesh(M.tri, x, y, u);
%title('Numerical solution of u(x,y)')
xlabel('x'); ylabel('y'); zlabel('u');
%xlim([-1 2]); ylim([-1 1]); view(2);

%Error of numerical solution
Uexact = [x y sin(pi * x) .* sin(pi * y)];
%Error = abs( u - Uexact(:,3) );
%Error_MAX = max(Error);

Error = norm(u - Uexact(:,3), Inf);

