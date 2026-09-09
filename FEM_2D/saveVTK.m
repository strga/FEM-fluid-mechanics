
function saveVTK(M, u, fname)

%Vrcholy
np = M.nbNod;

%Triangles
nt = M.nbTriangles;

% Save VTK datafor FE function u on mesh M
fid = fopen(fname, 'w+');
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Comment: Scalar data\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

% Vertices section
fprintf(fid, 'POINTS %d float\n', np);

for i = 1:np
    fprintf(fid, '%g %g %g\n', M.POS(i,1), M.POS(i,2), 0);
end

% CELLS WITH how many vertices + index list (VTK starts vertices with zero, thus shift -1 )
fprintf(fid, '\nCELLS %d %d\n', nt, 4 * nt);

% Intexy trojuhelniku
for i = 1:nt
    fprintf(fid, '3 %d %d %d\n', M.TRIANGLES(i, 1) - 1, M.TRIANGLES(i,2) - 1, M.TRIANGLES(i, 3) -1);
end

% CELL TYPE 5 - corresponding to VTK_TRIANGLE
fprintf(fid, '\n\nCELL_TYPES %d\n', nt);

for i = 1 : nt    
    fprintf(fid, '%d\n', 5);
end

%POINT DATA
fprintf(fid, '\nPOINT_DATA %d\n', np);
fprintf(fid, 'SCALARS %s float 1\n','u');
fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '%g\n', u);

%Uzavreni zapisu
fclose(fid);
end