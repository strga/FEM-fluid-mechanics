function s = computeIntegral3(M, f) %Středy strany

    s = 0;
    
    for k = 1:M.nbTriangles
        
        K = getMeshElement(M, k);

        Y = (K.B + K.C) / 2; Sa = f(Y(1), Y(2));
        Y = (K.C + K.A) / 2; Sb = f(Y(1), Y(2));
        Y = (K.A + K.B) / 2; Sc = f(Y(1), Y(2));

        str = (Sa+Sb+Sc)/3; %Strany
        
        s = s + K.vol * str;
    end
end