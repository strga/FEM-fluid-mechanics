
function [A, b] = Matice_tuhosti(M,f,Dir1,psi,alfa,ur)

I = [];
J = [];
VAL = [];
n = M.nbNod;
b = zeros(n, 1);

for k = 1:M.nbTriangles
    
    K = getMeshElement(M,k);
    G = computeGradient(K);
    ilist = [ K.iA, K.iB, K.iC ];
    
    Bk = eye(3,3);                                                   % Basis functions for coordinates of vertices of ref. triangle X1(0,0), X2(1,0), X3(0,1)...PHI_STRISKA_KSI,
    fk = [f(K.A(1), K.A(2)), f(K.B(1),K.B(2)), f(K.C(1),K.C(2))];    % Function evaluated in the vertices of triangle
    Omega = [1/3, 1/3, 1/3];                                         % Omega == vaha, |K| = h^2 / 2, |K| * Teziste, |K| * 1/3 * ( f(A) + f(B) + f(C) ), |K| * 1/3 * ( f(Sa) + f(Sb) + f(Sc) ),...
                                                                                ... Pro tesiste 1x == 1/2, pro vrcholy + stredy stran 3x == 1/6

    for i = 1:3
        for j = 1:3
            
            ak = K.vol * dot(G(:,i), G(:,j));
            I = [I, ilist(i)];
            J = [J, ilist(j)];
            VAL = [VAL, ak];  
            
        end
        
        b(M.tri(k,i)) = b(M.tri(k,i)) + K.vol * sum(Omega .* Bk(i,:) .* fk);
        
    end
end

A  = sparse(I, J, VAL); %Singularni matice

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%BOUNDARY CONDITIONS
Bmark = BoundaryPts(M); %Hranicni body

%NEUMANN BOUDARY CONDITION
iNEU = find(M.LINES(:, 3) == 10); %Find lines with given BC

for i = iNEU'
    S = getMeshBoundarySide(M, i);
    
    b(S.idxA) = b(S.idxA) + S.ds * psi( S.A(1), S.A(2), S.nn(1), S.nn(2) ) / 2; 
    b(S.idxB) = b(S.idxB) + S.ds * psi( S.B(1), S.B(2), S.nn(1), S.nn(2) ) / 2;
    
end
    
%ROBIN BOUNDARY CONDITION
iROB = find(M.LINES(:, 3) == 20); %Find lines with given BC

I = [];
J = [];
VAL = [];

for i = iROB'
    S = getMeshBoundarySide(M, i);
    
    b(S.idxA) = b(S.idxA) + S.ds * alfa * ur( S.A(1), S.A(2), S.nn(1), S.nn(2) ) / 2;
    b(S.idxB) = b(S.idxB) + S.ds * alfa * ur( S.B(1), S.B(2), S.nn(1), S.nn(2) ) / 2;
    
    I = [I, S.idxA, S.idxB];
    J = [J, S.idxA, S.idxB];
    
    hodnota = S.ds * alfa * 1 / 2;
    VAL = [VAL, hodnota, hodnota];
end

A = A + sparse(I, J, VAL, n, n);

%DIRICHLET BOUDARY CONDITION
iblist = find((Bmark == 1));

Pos1 = M.POS(iblist, 1:2); %Boundary points coordinates

b(iblist) = Dir1( Pos1(:,1), Pos1(:,2) );

%Internal points of domain
intlist = setdiff(1:n,iblist);

%Stiffness matrix - regular
A(iblist,:) = 0;
b(intlist) = b(intlist) - A(intlist,iblist) * b(iblist);
A(:,iblist) = 0;
A = A + sparse(iblist,iblist, 0 * iblist + 1, n, n );
end
