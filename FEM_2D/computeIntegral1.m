function s = computeIntegral1(M, f) %Těžiště

    s = 0;

    for k = 1:M.nbTriangles

        K = getMeshElement(M, k);

        T = (K.A + K.B + K.C)/3;

        s = s + K.vol*f(T(1),T(2));
    end
end
