function K = getMeshElement(M,k)        % Informace o konkrétním elementu

%Řádky 6...9 natahají do paměti struktury K informace z M.TRIANGLES ...
...- informace z k-tého řádku, kde na každém řádku jsou informace o ocislovanych vrcholech daneho trojuhelniku
        
K.iA = M.TRIANGLES(k,1);
K.iB = M.TRIANGLES(k,2);
K.iC = M.TRIANGLES(k,3);
K.mark = M.TRIANGLES(k,4);

K.A = M.POS(K.iA, 1:2)';     %Řádky 11...13 natahají do paměti informace o souradicich bodu trojuhelniku definovanem body K.iA...K.iC
K.B = M.POS(K.iB, 1:2)';
K.C = M.POS(K.iC, 1:2)';

BB = [K.B - K.A, K.C - K.A];    %BB je transformacni matice
K.BB = BB;                      %Přesunutí BB do struktury K jako K.BB, nic jiného se nezmění

K.inBB = inv(K.BB);             %Inverze transformační matice BB na inBB a přesonutá do struktury K.inBB

K.detBB = det(K.BB);

K.vol = det(K.detBB) / 2;            %Determinat je přesunut do struktury K, jako K.vol, vypocet plochy trojuhelnika

end
