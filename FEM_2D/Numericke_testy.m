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

x = M.POS(:,1); y = M.POS(:, 2);  %Souřadnice vrcholů elementů z M.POS
M.tri = M.TRIANGLES(1:M.nbTriangles,1:3); %Vytáhnout informace o trojuhelnicich

% figure(1)
% trisurf(M.tri, x, y, 2 .* pi.^2 .* sin(pi*x) .* sin(pi*y));
% title('Function f(x,y)') 

%Prava strana Poissonovy rovnice
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

%{
%CHYBY PRO SQUARE_DIR
DIR.h4 = 0.053029287545609; ...5x5
DIR.h8 = 0.012950746721908; ...9x9
DIR.h16 = 0.003218964440088; ...17x17
DIR.h32 = 0.0008035776793775540; ...33x33

%CHYBA NA OBLASTI SQUARE_DIR_NEU1
NEU.h4 = 0.049528930575971; ...5x5
NEU.h8 = 0.012676066895710; ...9x9
NEU.h16 = 0.003192452264879; ...17x17
NEU.h32 = 0.0007996691096511714; ...33x33


%CHYBA NA OBLASTI SQUARE_DIR_ROB1
ROB.h4 = 0.044921387047463; ...5x5
ROB.h8 = 0.010991997220743; ...9x9
ROB.h16 = 0.002766567045471; ...17x17
ROB.h32 = 0.0006909597353929042; ...33x33

%CHYBA NA OBLASTI SQUARE_DIR_NEU_ROB
MIX.h4 = 0.080924522037904; ...5x5
MIX.h8 = 0.021272515399067; ...9x9
MIX.h16 = 0.005374663442235; ...17x17
MIX.h32 = 0.001348237778049; ...33x33

DIR.E21 = ((DIR.h8) / (DIR.h4)).^(1);
DIR.E32 = ((DIR.h16) / (DIR.h8)).^(1);
DIR.E43 = ((DIR.h32) / (DIR.h16)).^(1);

DIR.p21 = -log2(DIR.E21);
DIR.p32 = -log2(DIR.E32);
DIR.p43 = -log2(DIR.E43);

NEU.E21 = ((NEU.h8) / (NEU.h4)).^(1);
NEU.E32 = ((NEU.h16) / (NEU.h8)).^(1);
NEU.E43 = ((NEU.h32) / (NEU.h16)).^(1);

NEU.p21 = -log2(NEU.E21);
NEU.p32 = -log2(NEU.E32);
NEU.p43 = -log2(NEU.E43);

ROB.E21 = ((ROB.h8) / (ROB.h4)).^(1);
ROB.E32 = ((ROB.h16) / (ROB.h8)).^(1);
ROB.E43 = ((ROB.h32) / (ROB.h16)).^(1);

ROB.p21 = -log2(ROB.E21);
ROB.p32 = -log2(ROB.E32);
ROB.p43 = -log2(ROB.E43);

MIX.E21 = ((MIX.h8) / (MIX.h4)).^(1);
MIX.E32 = ((MIX.h16) / (MIX.h8)).^(1);
MIX.E43 = ((MIX.h32) / (MIX.h16)).^(1);

MIX.p21 = -log2(MIX.E21);
MIX.p32 = -log2(MIX.E32);
MIX.p43 = -log2(MIX.E43);

%Chyba Semicircle
%Dirichlet
    %triangles = 2017
    %Emax = 0.010525774606005
    %BoundaryPoints =  1085

%Neumann =
    %triangles = 2017
    %Emax = 0.017052166058456
    %BoundaryPoints = 1085
    
%Robin = 
    %triangles = 2017
    %Emax = 0.018263584031156
    %BoundaryPoints = 1085
%}