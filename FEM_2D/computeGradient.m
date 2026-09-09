function RR = computeGradient(K)

GH = [-1 -1; 1 0; 0 1]; % Derivace bazove funkce na referencnim elementu

RR = (GH * K.inBB)'; % Matrix 2x3, 

end