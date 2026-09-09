function K = getMeshElement(M,k)        % Informace o konkrétním elementu

%Vertices and marks of triangle elements
        
K.iA = M.TRIANGLES(k,1);
K.iB = M.TRIANGLES(k,2);
K.iC = M.TRIANGLES(k,3);
K.mark = M.TRIANGLES(k,4);

K.A = M.POS(K.iA, 1:2)';     % Coordinates of vertices
K.B = M.POS(K.iB, 1:2)';
K.C = M.POS(K.iC, 1:2)';

BB = [K.B - K.A, K.C - K.A];    %BB - transformation matrix
K.BB = BB;

K.inBB = inv(K.BB);             %Inverze transfor of matrix

K.detBB = det(K.BB);

K.vol = det(K.detBB) / 2;            %Area of element / triangle

end
