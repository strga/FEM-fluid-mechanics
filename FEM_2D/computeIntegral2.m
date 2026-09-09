function s = computeIntegral2(M, f) %Vrcholy

    s = 0;

    for k = 1:M.nbTriangles
    
        K = getMeshElement(M, k);

        a = f(K.A(1),K.A(2));
        b = f(K.B(1),K.B(2));
        c = f(K.C(1),K.C(2));

        V = (a+b+c)/3; %Vrcholy
        
        s = s + K.vol * V;
    end
end