clear all
close all
format long

M = load_gmsh('NACA_0012.msh');
%M = load_gmsh('Zukovskij.msh');

x = M.POS(:,1); y = M.POS(:, 2);
M.tri = M.TRIANGLES(1:M.nbTriangles,1:3);

theta = 0;  %[deg]
u.U = 50;     %[m/s]
rho = 1.225;  %Hustota vzduchu [kg/m^3]
max_speed = 200;

u1 = u.U * cosd(theta);
u2 = u.U * sind(theta);

%Right side of Poisson
bvp.f = @(x, y) 0 * x;

%Dirchlet boundary condition - wing
bvp.MarkDir1 = 10;
bvp.Dir1 = @(x,y) -0.030811632588222;

%Dirchlet boundary condition - input
bvp.MarkDir2 = 2;
bvp.Dir2 = @(x, y) (u1 * y - u2 * x);

%Neumann boundary condition - output
bvp.MarkNeu = 10;
bvp.phiNeu = @(x, y, n1, n2) ( n2 * u1 - n1 * u2 );

%INCOMPRESSIBLE FLOW
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% [A,b] = DiscretizePotentialFlow(M,bvp);
% 
% %Potential flow
% psi = A \ b;

% COMPRESSIBLE FLOW
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
[density.U1, density.rho_grad] = density_gradientPSI(max_speed, rho);
figure(1)
plot(density.U1, density.rho_grad); %Relationship speed [grad(U)] and density
xlabel('Rychlost'); ylabel('Hustota');
% Initial computation
u1 = zeros(M.nbNod,1);
[A,b] = DiscretizePotentialFlowCompressible(M, u1, bvp, density);
psi = A \ b;
Error = norm(psi - u1); %Initial error
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

%Find value of Stream function (Proudova fce)
while (Error > 1e-9)
   
    psi_in = psi;
    
    [A,b] = DiscretizePotentialFlowCompressible(M, psi_in, bvp, density);
    psi = A \ b;
    
    Error = norm(psi - psi_in);

end

% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

[MM,bbx,bby] = Matice_rekonstrukce(M,psi);

%Speed vector elements
u.uux = MM \ bby;
u.uuy = -(MM \ bbx);

[V,p,Cp,Profile] = physical_quantities(M,rho,u);

Edge = find( M.LINES(:,3) == 1 );

FD = 0;
FL = 0;

for i = Edge'
    
    S = getMeshBoundarySide(M, i);
    
    p_avg = mean( [ p.wing(S.idxA), p.wing(S.idxB) ] );
    
    FD = FD + p_avg * S.ds * S.nn(1); %Drag force
    
    FL = FL + p_avg * S.ds * S.nn(2); %Lift force
end

% For Cd and Cl the wing has unit area
Cd = FD / ( 0.5 * rho * u.U^2); %Drag
Cl = FL / ( 0.5 * rho * u.U^2); %Lift
% Approximation of lift coefficient Cl
CL = 2 * pi * sind( theta );

%Finding value of constant (Wing Dirichlet B.C.)
[Constants,Psi_point] = Plot_constant(M,bvp);

%Extrapolation of CONSTANT for Dir. B.C.
CONSTANT = interp1(Psi_point - Constants',Constants,0,'linear','extrap');
fprintf('Dirichlet constant = %.5cf\n', CONSTANT);

figure(2)
plot(Constants,Psi_point - Constants') %Value of constant is the x value of the intercept of the line and yline(0)
hold on
yline(0)
hold off

%Sousedni elementy k uzlu na konci kridla
End.vertex = find( M.POS(:,1) == max(Profile.length) );
End.vert_coord = M.POS(End.vertex, 1:2);
End.elements = find( M.TRIANGLES(:,1) == End.vertex | M.TRIANGLES(:,2) == End.vertex | M.TRIANGLES(:,3) == End.vertex );
End.points = unique( M.TRIANGLES( End.elements, 1:3 ) );
End.coord = M.POS( End.points, 1:2 );

%{
figure(3) %Proudova fce
trimesh(M.tri, x, y, psi);
xlabel('x'); ylabel('y'); zlabel('u');
title('Potential flow')
xlim([-1 2]); ylim([-1 1]); view(2);

figure(4) %Vektorove pole rychlosti
quiver(x,y,u.uux,u.uuy,0.5);
title('Vector field')
xlim([-1 2]); ylim([-1 1]); view(2);
%}

%VEKTOROVE POLE --- pro Zukovskeho profil - kontrola
figure(99)
quiver(x(980:1000),y(980:1000),u.uux(980:1000),u.uuy(980:1000))

%VEKTOROVE POLE --- KONTROLA SPLNENI KUTTA-ZUK. PODMINKY
figure(100)
quiver(End.coord(:,1), End.coord(:,2), u.uux(End.points), u.uuy(End.points), 0.5);
title('Vector field')
hold on
plot(1,0,'r*'), plot(1.001,0,'r*');
hold off

%{
%Pressure around airfoil
figure(5)
plot(Profile.length, p.upper)
title('Pressure distribution upper part')
xlabel('Length of wing'); ylabel('Pressure');

figure(6)
plot(Profile.length, p.lower)
title('Pressure distribution lower part')
xlabel('Length of wing'); ylabel('Pressure');
%}

figure(7)
plot(Profile.length, p.rozdil)
title('Pressure difference')
xlabel('Délka křídla'); ylabel('Rozdíl tlaku');

%{
%Speed aroud airfoil
figure(8)
plot(Profile.length, V.upper)
title('Speed along upper part')
xlabel('Length of wing'); ylabel('Speed');

figure(9)
plot(Profile.length, V.lower)
title('Speed along lower part')
xlabel('Length of wing'); ylabel('Speed');
%}

figure(10)
plot(Profile.length, Cp.lower,Profile.length,Cp.upper)
%title('Cp')
xlabel('Délka křídla'); ylabel('Cp');
legend({'Spodní povrch','Horní povrch'},'Location','southwest')

figure(11)
plot(Profile.length, Cp.diff)
%title('Cp')
xlabel('Délka křídla'); ylabel('Rozdíl Cp');
%legend({'Spodní povrch','Horní povrch'},'Location','southwest')

%SPEED IN AREA OMEGA

%{
figure(12)
subplot(1,2,1);
h = trisurf(M.tri, x, y, u.uux);
set(h,'EdgeAlpha', 0);
title('Speed x');
xlim([-1 1]); ylim([-1 1]); view(2);
colorbar
axis equal
subplot(1,2,2);
%figure(8)
h = trisurf(M.tri, x, y, u.uuy);
set(h,'EdgeAlpha', 0);
title('Speed y');
xlim([-1 1]); ylim([-1 1]); view(2);
colorbar
axis equal
hold off
%}

figure(13) % Speed
h = trisurf(M.tri, x, y, V.Omega);
set(h,'EdgeAlpha', 0)
%title('Speed')
xlim([-0.25 1.25]); ylim([-0.5 0.5]); view(2);
colorbar

%PRESSURE IN OMEGA

figure(14) %Pressure
h = trisurf(M.tri, x, y, p.Omega);
set(h,'EdgeAlpha', 0);
%title('Pressure');
view(2);
%daspect([3 2 1])
xlim([-0.25 1.25]); ylim([-0.5 0.5]);
colorbar

%{
%Potential flow / proudova fce
figure(15)
trimesh(M.tri, x, y, psi)
xlabel('x'); ylabel('y'); zlabel('u');
title('Potential flow')
%xlim([0.98 1.02]); ylim([-0.02 0.02]);
view(2);
hold on

figure(16) %rychlost smer X
h = trisurf(M.tri, x, y, u.uux);
set(h,'EdgeAlpha', 0);
title('Speed x - detail');
xlim([0.98 1.02]); ylim([-0.02 0.02]); view(2);
colorbar

figure(17) %rychlosti smer Y
h = trisurf(M.tri, x, y, u.uuy);
set(h,'EdgeAlpha', 0);
title('Speed y - detail');
xlim([0.98 1.02]); ylim([-0.02 0.02]); view(2);
colorbar

figure(18) %rychlost velikost
h = trisurf(M.tri, x, y, V.Omega);
set(h,'EdgeAlpha', 0)
title('Speed VELIKOST - detail')
xlim([0.98 1.02]); ylim([-0.02 0.02]); view(2);
colorbar
%}

%Compress flow

% Rho.Omega = interp1(density.U1, density.rho_grad, V.Omega);
% 
% figure(19) %Density distribution compress. flow
% h = trisurf(M.tri, x, y, Rho.Omega);
% set(h, 'EdgeAlpha', 0);
% view(2)
% xlim([-0.25 1.25]); ylim([-0.5 0.5]);
% colorbar
