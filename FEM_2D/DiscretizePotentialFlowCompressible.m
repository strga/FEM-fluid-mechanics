function [A,b] = DiscretizePotentialFlowCompressible(M,u1,bvp,density)

I = [];
J = [];
VAL = [];
n = M.nbNod;
b = zeros(n, 1);

velPsi = density.U1; % Speed
rhoPsi = density.rho_grad; % Density for speed U1

for k = 1:M.nbTriangles
    
    K = getMeshElement(M,k);
    G = computeGradient(K);
    ilist = [ K.iA, K.iB, K.iC ];
    
    %Compressible addition
    grad_uK = u1(K.iA) * G(:,1) + u1(K.iB) * G(:,2) + u1(K.iC) * G(:,3); %Gradient PSI
    size_grad_uK = norm(grad_uK); %Size |Grad(PSI)|
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    Bk = eye(3,3);
    fk = [bvp.f(K.A(1), K.A(2)), bvp.f(K.B(1),K.B(2)), bvp.f(K.C(1),K.C(2))];   % funkcni hodnoty ve vrcholech trojuhelniku
    Omega = [1/3, 1/3, 1/3]; %Vaha
    
    %Inverse density at point size_grad_uK
    rhoK = interp1(velPsi, rhoPsi, size_grad_uK);
    %Inverse density
    invrhoK = 1 / rhoK;                                                                        
    
    
    for i = 1:3
        for j = 1:3
            
            ak = K.vol * invrhoK * dot(G(:,i), G(:,j)); %Addition of inverse density
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
Bmark = BoundaryPts(M);

%NEUMANN BOUDARY CONDITION
iNEU = find(M.LINES(:, 3) == 10); %Find lines with given BC

for i = iNEU'
    S = getMeshBoundarySide(M, i);
    
    b(S.idxA) = b(S.idxA) + S.ds * bvp.phiNeu( S.A(1), S.A(2), S.nn(1), S.nn(2) ) / 2; 
    b(S.idxB) = b(S.idxB) + S.ds * bvp.phiNeu( S.B(1), S.B(2), S.nn(1), S.nn(2) ) / 2;
    
end

%DIRICHLET BOUDARY CONDITION
iblist = find((Bmark == 1) | (Bmark == 2));
iblist1 = find((Bmark == 1));
iblist2 = find((Bmark == 2));

Pos1 = M.POS(iblist1, 1:2); %Boundary points coordinates
Pos2 = M.POS(iblist2, 1:2); %Boundary points coordinates

b(iblist1) = bvp.Dir1( Pos1(:,1), Pos1(:,2) );
b(iblist2) = bvp.Dir2( Pos2(:,1), Pos2(:,2) );

%Internal points of domain
intlist = setdiff(1:n,iblist);

%Matice tuhosti - regularni
A(iblist,:) = 0;
b(intlist) = b(intlist) - A(intlist,iblist) * b(iblist);
A(:,iblist) = 0;
A = A + sparse(iblist,iblist, 0 * iblist + 1, n, n );

end
