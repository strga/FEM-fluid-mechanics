function RR = computeGradient(K)

GH = [-1 -1; 1 0; 0 1]; % Derivatives of basis function on reference element

RR = (GH * K.inBB)'; % Matrix 2x3, 

end
