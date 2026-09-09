function s = computeIntegral1(M, f) %Těžiště

    s = 0;

    for k = 1:M.nbTriangles

        K = getMeshElement(M, k);

        %Sab = (K.A + K.B)/2;

        T = (K.A + K.B + K.C)/3;
        %T = Sab + 1/3 * (K.C - Sab);

        s = s + K.vol*f(T(1),T(2));
    end
end