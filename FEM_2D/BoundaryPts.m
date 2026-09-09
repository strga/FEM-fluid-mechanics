function Bmark = BoundaryPts(M)

Bmark = zeros(M.nbNod, 1);

%Find line with given Physical values

%Neummann - find lines on the border of region with given B.C.
k = find(M.LINES(:, 3) == 10);
Points3 = unique([M.LINES(k,1); M.LINES(k,2)]);

%Robin - find lines on the border of region with given B.C.
k = find(M.LINES(:, 3) == 20);
Points2 = unique([M.LINES(k,1); M.LINES(k,2)]);

%Dirichlet - find lines on the border of region with given B.C.
k = find(M.LINES(:, 3) == 1);
Points1 = unique([M.LINES(k,1); M.LINES(k,2)]);

%Dirichlet - find lines on the border of region with given B.C.
k = find(M.LINES(:, 3) == 2);
Points0 = unique([M.LINES(k,1); M.LINES(k,2)]);

%Neumann == 10
Bmark(Points3) = 10;

%Robin == 20
Bmark(Points2) = 20;

%Dirichlet == 1
Bmark(Points1) = 1;

%Dirichlet == 2
Bmark(Points0) = 2;

end