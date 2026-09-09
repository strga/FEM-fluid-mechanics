function vol = getDomainVolume(M)       % Výpočet plochy

    vol = 0;
    
    for k = 1:M.nbTriangles
        K = getMeshElement(M, k);
        vol = vol + K.vol;
    end
end
