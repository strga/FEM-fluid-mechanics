
function [MM, bbx, bby] = Matice_rekonstrukce(M,psi)

I = [];
J = [];
VAL = [];
n = M.nbNod;
bbx = zeros(n, 1);
bby = zeros(n, 1);

% Stred stran trojuhelniku S
W = [0 1/2 1/2; 
    1/2 0 1/2; 
    1/2 1/2 0];

for k = 1:M.nbTriangles
    
    K = getMeshElement(M,k);
    G = computeGradient(K);
    psi_i = [psi(K.iA), psi(K.iB), psi(K.iC)]; %Instead of index list (ilist)
    
    for i = 1:3
        for j = 1:3
   
            ak = 1/3 * K.vol * dot( W(i,:), W(j,:) ); %
            
            I = [I; M.tri(k,i)];
            J = [J; M.tri(k,j)];
            VAL = [VAL; ak];
            
        end
        
        %Partial derivetive
        Psi_der_x = dot( psi_i, G(1,:) );
        Psi_der_y = dot( psi_i, G(2,:) );
        
        bbx(M.tri(k,i)) = bbx(M.tri(k,i)) + 1/3 * K.vol * Psi_der_x;
        bby(M.tri(k,i)) = bby(M.tri(k,i)) + 1/3 * K.vol * Psi_der_y;
        
    end
end

MM  = sparse(I, J, VAL);
%{
%BOUNDARY CONDITIONS
Bmark = BoundaryPts(M); %Hranicni body

%NEUMANN BOUDARY CONDITION
iNEU = find(M.LINES(:, 3) == 10); %Find lines with given BC

for i = iNEU'
    S = getMeshBoundarySide(M, i);
    
    b(S.idxA) = b(S.idxA) + S.ds * psi( S.A(1), S.A(2), S.nn(1), S.nn(2) ) / 2; 
    b(S.idxB) = b(S.idxB) + S.ds * psi( S.B(1), S.B(2), S.nn(1), S.nn(2) ) / 2;
    
end

%DIRICHLET BOUDARY CONDITION
%iblist = find((Bmark == 1));
iblist = find((Bmark == 1) | (Bmark == 2));
iblist1 = find((Bmark == 1));
iblist2 = find((Bmark == 2));

Pos1 = M.POS(iblist1, 1:2); %Boundary points coordinates
Pos2 = M.POS(iblist2, 1:2); %Boundary points coordinates

b(iblist1) = Dir1( Pos1(:,1), Pos1(:,2) );
b(iblist2) = Dir2( Pos2(:,1), Pos2(:,2) );

%Internal points of domain
intlist = setdiff(1:n,iblist);

%Matice tuhosti - regularni
A(iblist,:) = 0;
b(intlist) = b(intlist) - A(intlist,iblist) * b(iblist);
A(:,iblist) = 0;
A = A + sparse(iblist,iblist, 0 * iblist + 1, n, n );
%}
end